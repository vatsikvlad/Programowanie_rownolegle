#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 16
#define G6_LEN 21 
#define BATCH_SIZE 500000 

__global__ void check_integral_graphs(const char* d_strings, int* d_results, int num_graphs) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_graphs) return;

    d_results[tid] = 0; 

    // odczyt tekstu
    const char* my_str = d_strings + (tid * G6_LEN);
    
    if (my_str[0] != 'O') return;

    double a[N*(N-1)/2 + N + 1];
    
    // dekodwanie
    int p = 1; 
    int cur_bit = 0; 
    unsigned char cur_byte = 0;

    for (int j = 0; j < N; j++) {
        for (int i = 0; i < j; i++) {
            if (cur_bit == 0) {
                cur_byte = my_str[p++] - 63;
                cur_bit = 6;
            }
            int edge = (cur_byte >> (cur_bit - 1)) & 1;
            cur_bit--;

            int idx = (j * (j + 1) / 2) + i + 1;
            a[idx] = (double)edge; 
        }
        int diag_idx = (j * (j + 1) / 2) + j + 1;
        a[diag_idx] = 0.0;
    }

    // householder
    int i, j, k, k3, k4, L, L1, z, cond;
    double eps, g, h, ma, mn, norm, s, t, u, w;
    double d[N+1], e[N+1], e2[N+1], Lb[N+1], x[N+1];

    i = 0;
    for (L=1;L<=N;L++) { i += L; d[L] = a[i]; }

    for (L=N;L>=2;L--) {
        i--; j = i; h = a[j]; s = 0;
        for (k=L-2;k>=1;k--) { i--; g = a[i]; s += g*g; }
        i--;
        if (s == 0) { e[L] = h; e2[L] = h*h; a[j] = 0.0; }
        else {
            s += h*h; e2[L] = s; g = sqrt(s); if (h>=0.0) g=-g;  
            e[L] = g;
            s = 1.0 / (s-h*g);
            a[j] = h - g; h = 0.0; L1 = L - 1; k3 = 1;
            for (j=1;j<=L1;j++) {
                k4 = k3; g = 0;
                for (k=1;k<=L1;k++) { 
                    g +=a[k4]*a[i+k]; 
                    if (k<j)  z = 1; else z = k;
                    k4 += z; 
                }
                k3 += j; g *= s; e[j] = g; h += a[i+j]*g;
            }
            h *= 0.5*s; k3 = 1;
            for (j=1;j<=L1;j++) {
                s = a[i+j]; g = e[j]-h*s; e[j] = g;
                for (k=1;k<=j;k++) { a[k3] += -s*e[k]-a[i+k]*g; k3++; }
            }
        }
        h = d[L]; d[L] = a[i+L]; a[i+L] = h;
    }
    h = d[1]; d[1] = a[1]; a[1] = h; e[1] = 0.0; e2[1] = 0.0; s = d[N];
    t = fabs(e[N]); mn = s - t; ma = s + t; 
    
    for (i=N-1;i>=1;i--) {
        u = fabs(e[i]); h = t + u; t = u; s = d[i]; u = s - h;
        if (u < mn) mn = u;
        u = s + h;
        if (u > ma) ma = u;
    }
    for (i=1;i<=N;i++) { Lb[i] = mn; x[i] = ma; }
    norm = fabs(mn); s = fabs(ma);
    if (s>norm) norm = s;
    w = ma; eps = 7.28e-17*norm;

    for (k=N;k>=1;k--) {
        s = mn; i = k;
        do {
            cond = 0; g = Lb[i];
            if (s < g) s = g; else { i--; if (i>=1) cond = 1; }
        } while (cond);
        g = x[k];
        if (w>g) w = g;
        
        while (w-s>2.91e-16*(fabs(s)+fabs(w))+eps) { 
            if (floor(w+1e-5)<s-1e-5) { d_results[tid] = 0; return; } 
            L1 = 0; g = 1.0; t = 0.5*(s+w);
            for (i=1;i<=N;i++) {
                if (g!=0)  g = e2[i] / g; else g = fabs(6.87e15*e[i]); 
                g = d[i] - t - g; 
                if (g<0) L1++;
            }
            if (L1<1) { s = t; Lb[1] = s; }
            else { 
                if (L1<k) {
                    s = t; Lb[L1+1] = s;
                    if (x[L1]>t) x[L1] = t;
                }
                else w = t;
            }
        }
        u = 0.5*(s+w); x[k] = u;
        if  (!(( ceil(u) - u  < 1e-5 ) || ( u - floor(u) < 1e-5 ))) { d_results[tid] = 0; return; } 
    }
    
    d_results[tid] = 1; 
}


int main(int argc, char* argv[]) {
    int block_size = 256; 
    
    char* h_strings = (char*)malloc(BATCH_SIZE * G6_LEN);
    int* h_results = (int*)malloc(BATCH_SIZE * sizeof(int));

    char* d_strings;
    int* d_results;
    cudaMalloc((void**)&d_strings, BATCH_SIZE * G6_LEN);
    cudaMalloc((void**)&d_results, BATCH_SIZE * sizeof(int));

    int current_batch_size = 0;
    char buffer[1024]; 

    while (fgets(buffer, sizeof(buffer), stdin)) {
        if (buffer[0] != 'O') continue; 
        
        memcpy(h_strings + (current_batch_size * G6_LEN), buffer, G6_LEN);
        current_batch_size++;

        if (current_batch_size == BATCH_SIZE) {
            
            cudaMemcpy(d_strings, h_strings, BATCH_SIZE * G6_LEN, cudaMemcpyHostToDevice);
            
            int num_blocks = (BATCH_SIZE + block_size - 1) / block_size;
            check_integral_graphs<<<num_blocks, block_size>>>(d_strings, d_results, BATCH_SIZE);
            cudaDeviceSynchronize();
            
            cudaMemcpy(h_results, d_results, BATCH_SIZE * sizeof(int), cudaMemcpyDeviceToHost);
            
            for (int i = 0; i < BATCH_SIZE; i++) {
                if (h_results[i] == 1) {
                    char out_buf[G6_LEN + 2];
                    memcpy(out_buf, h_strings + (i * G6_LEN), G6_LEN);
                    out_buf[G6_LEN] = '\n';
                    out_buf[G6_LEN + 1] = '\0';
                    fputs(out_buf, stdout);
                }
            }
            current_batch_size = 0;
        }
    }

    // sprzątanie
    if (current_batch_size > 0) {
        cudaMemcpy(d_strings, h_strings, current_batch_size * G6_LEN, cudaMemcpyHostToDevice);
        int num_blocks = (current_batch_size + block_size - 1) / block_size;
        check_integral_graphs<<<num_blocks, block_size>>>(d_strings, d_results, current_batch_size);
        cudaDeviceSynchronize();
        cudaMemcpy(h_results, d_results, current_batch_size * sizeof(int), cudaMemcpyDeviceToHost);
        
        for (int i = 0; i < current_batch_size; i++) {
            if (h_results[i] == 1) {
                char out_buf[G6_LEN + 2];
                memcpy(out_buf, h_strings + (i * G6_LEN), G6_LEN);
                out_buf[G6_LEN] = '\n';
                out_buf[G6_LEN + 1] = '\0';
                fputs(out_buf, stdout);
            }
        }
    }

    free(h_strings);
    free(h_results);
    cudaFree(d_strings);
    cudaFree(d_results);

    return EXIT_SUCCESS;
}