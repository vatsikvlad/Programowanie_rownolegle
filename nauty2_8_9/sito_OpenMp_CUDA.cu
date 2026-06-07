/** 
@file=sito8.cu
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
 
#define NMAX 20
#define GLEN 21
#define BATCH_SIZE 65536
 
__global__ void test(char *d_bus, int *d_results, int total_graphs) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total_graphs) return;

    char buffer[GLEN];
    for(int i = 0; i < GLEN; i++) {
        buffer[i] = d_bus[idx * GLEN + i];
    }
 
    int i, j, k, k3, k4, L, L1, z;
    double eps, g, h, ma, mn, norm, s, t, u, w;
    int cond;
    double d[NMAX + 1], e[NMAX + 1], e2[NMAX + 1], Lb[NMAX + 1];
    double x[NMAX + 1];
    double a[NMAX * (NMAX - 1) / 2 + NMAX + 1];
    int n;
    int bit, poz, poz2;
    bit = 32;
    poz = 1;
    poz2 = 1;

    char *BUFFOR = buffer;

    n = BUFFOR[0] - 63;
    a[0] = 0.0;
    for (i = 0; i < n; i++){
        for (j = 0; j <= i; j++)
        {
            if (i == j) {
                a[poz2++] = 0; 
            }
            else {
                if (bit == 0) { 
                    bit = 32;  
                    poz++; 
                }
                if ((BUFFOR[poz] - 63) & bit){ 
                    a[poz2++] = 1; 
                }
                else{ 
                    a[poz2++] = 0; 
                }
                bit = bit >> 1;
            }
        }
    }
    
    int k1 = 1;
    int k2 = n;
    i = 0;
    for (L = 1; L <= n; L++) { 
        i += L; 
        d[L] = a[i]; 
    }

    for (L = n; L >= 2; L--)
    {
        i--; 
        j = i; 
        h = a[j]; 
        s = 0;
        for (k = L - 2; k >= 1; k--) { 
            i--; 
            g = a[i]; 
            s += g * g; 
        }
        i--;
        if (s == 0) {
            e[L] = h; 
            e2[L] = h * h; 
            a[j] = 0.0; 
        }
        else {
            s += h * h; 
            e2[L] = s; 
            g = sqrtf(s); 
            if (h >= 0.0) {
                g=-g;
            }
            e[L] = g;
            s = 1.0 / (s - h * g);
            a[j] = h - g; 
            h = 0.0; 
            L1 = L - 1; 
            k3 = 1;
            for (j = 1; j <= L1; j++)
            {
                k4 = k3; g = 0;
                for (k = 1; k <= L1; k++) { 
                    g += a[k4] * a[i + k]; 
                    if (k < j){
                        z = 1;
                    }
                    else{
                        z = k;
                    }
                    k4 += z; 
                }
                k3 += j; 
                g *= s; 
                e[j] = g; 
                h += a[i + j] * g;
            }
            h *= 0.5 * s; 
            k3 = 1;
            for (j = 1; j <= L1; j++)
            {
                s = a[i + j]; 
                g = e[j] - h * s; 
                e[j] = g;
                for (k = 1; k <= j; k++) { 
                    a[k3] += -s * e[k] - a[i + k] * g; 
                    k3++; 
                }
            }
        }
        h = d[L]; 
        d[L] = a[i + L]; 
        a[i + L] = h;
    }
    h = d[1]; 
    d[1] = a[1]; 
    a[1] = h; 
    e[1] = 0.0; 
    e2[1] = 0.0; 
    s = d[n];
    t = fabsf(e[n]); 
    mn = s - t; 
    ma = s + t;
    for (i = n - 1; i >= 1; i--)
    {
        u = fabsf(e[i]); 
        h = t + u; 
        t = u; 
        s = d[i]; 
        u = s - h;
        if (u < mn) {
            mn = u;
        }
        u = s + h;
        if (u > ma) {
            ma = u;
        }
    }
    for (i = 1; i <= n; i++) { 
        Lb[i] = mn; 
        x[i] = ma; 
    }
    norm = fabsf(mn); 
    s = fabsf(ma);
    if (s > norm) {
        norm = s;
    }
    w = ma; 
    eps = 7.28e-17f * norm;
    for (k = k2; k >= k1; k--)
    {
        s = mn; 
        i = k;
        do {
            cond = 0; 
            g = Lb[i];
            if (s < g) {
                s = g;
            } 
            else { 
                i--; 
                if (i >= k1) {
                    cond = 1;
                } 
            }
        } while (cond);
        g = x[k];
        if (w > g) {
            w = g;
        }
        while (w - s > 2.91e-16 * (fabsf(s) + fabsf(w)) + eps)
        {
            if (floorf(w + 10e-5) < s - 10e-5) {
                return;
            }  // przedział nie zawiera liczby całkowitej
            L1 = 0; 
            g = 1.0; 
            t = 0.5 * (s + w);
            for (i = 1; i <= n; i++)
            {
                if (g != 0) {
                    g = e2[i] / g;
                } 
                else {
                    g = fabsf(6.87e15 * e[i]);
                }
                g = d[i] - t - g;
                if (g < 0) {
                    L1++;
                }
            }
            if (L1 < k1) { 
                s = t; 
                Lb[k1] = s; 
            }
            else
            { 
                if (L1 < k)
                {
                    s = t; 
                    Lb[L1 + 1] = s;
                    if (x[L1] > t) {
                        x[L1] = t;
                    }
                }
                else {
                    w = t;
                }
            }
        }
        u = 0.5 * (s + w); 
        x[k] = u;
        if (!((ceilf(u) - u < 10e-5) || (u - floorf(u) < 10e-5))) { 
            return; 
        }
    }
    d_results[idx] = 1;
}

int main(int argc, char *argv[])
{
    char *h_batch = (char*)malloc(BATCH_SIZE * GLEN);
    int *h_results = (int*)malloc(BATCH_SIZE * sizeof(int));

    char *d_batch;
    int *d_results;
    cudaMalloc((void**)&d_batch, BATCH_SIZE * GLEN);
    cudaMalloc((void**)&d_results, BATCH_SIZE * sizeof(int));

    char line[1024];
    int count = 0;

    while (fgets(line, sizeof(line), stdin)) {
        __builtin_memcpy(&h_batch[count * GLEN], line, GLEN);
        h_results[count] = 0;
        count++;       

        if (count == BATCH_SIZE) {
            cudaMemcpy(d_batch, h_batch, BATCH_SIZE * GLEN, cudaMemcpyHostToDevice);
            cudaMemcpy(d_results, h_results, BATCH_SIZE * sizeof(int), cudaMemcpyHostToDevice);

            int threadsPerBlock = 256;
            int blocksPerGrid = (BATCH_SIZE + threadsPerBlock - 1) / threadsPerBlock;

            test<<<blocksPerGrid, threadsPerBlock>>>(d_batch, d_results, BATCH_SIZE);
            
            cudaMemcpy(h_results, d_results, BATCH_SIZE * sizeof(int), cudaMemcpyDeviceToHost);

            #pragma omp parallel for schedule(static) // robimy równoległy cykl i dzielimy iteracje na równe bloki między wątkami
            for (int i = 0; i < BATCH_SIZE; i++) {
                if (h_results[i] == 1) { //jeżeli znalezliśmy graf to wchodzimy do if
                    char out_buf[GLEN + 2]; //tworzymy tymczasowy bufor dla wątku
                    __builtin_memcpy(out_buf, &h_batch[i * GLEN], GLEN); //kopijujemy graf do lokalnego buforu
                    out_buf[GLEN] = '\n'; // dodajemy na końcu konic linii
                    out_buf[GLEN + 1] = '\0'; // dodajemy zero aby fputs wiedział, gdzie się kończy wiersz w pamięci

                    #pragma omp critical // Tylko jeden wątek może wywołać poniższą funkcję
                    {
                        fputs(out_buf, stdout); // Wypisujemy znaleziony graf
                    }
                }
            }
            count = 0;
        }
    }

    if (count > 0) {
        cudaMemcpy(d_batch, h_batch, count * GLEN, cudaMemcpyHostToDevice);
        cudaMemcpy(d_results, h_results, count * sizeof(int), cudaMemcpyHostToDevice);
        int threadsPerBlock = 256;
        int blocksPerGrid = (count + threadsPerBlock - 1) / threadsPerBlock;
        
        test<<<blocksPerGrid, threadsPerBlock>>>(d_batch, d_results, count);
        cudaDeviceSynchronize();
        cudaMemcpy(h_results, d_results, count * sizeof(int), cudaMemcpyDeviceToHost);

        #pragma omp parallel for schedule(static) //wykonujemy tutaj to samo co i w górze, ale dla pozostałych grafów w pakiecie
        for (int i = 0; i < count; i++) {
            if (h_results[i] == 1) {
                char out_buf[GLEN + 2];
                __builtin_memcpy(out_buf, &h_batch[i * GLEN], GLEN);
                out_buf[GLEN] = '\n';
                out_buf[GLEN + 1] = '\0';
                
                #pragma omp critical
                {
                    fputs(out_buf, stdout);
                }
            }
        }
    }

    cudaFree(d_batch);
    cudaFree(d_results);
    free(h_batch);
    free(h_results);

    return EXIT_SUCCESS;
}