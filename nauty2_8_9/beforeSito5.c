/**
@file sito5.c 
 
Kompilacja: 
gcc -O3 sito5.c -o sito5 -lm
 
$ gcc -O3 sito5.c -o sito5 -fopenmp -lm
 
Uruchomienie (zliczanie bez zapamiętywania znalezionych grafów):
 
student@ANT:~/nauty/nauty27r3$ ./geng -c 5 2>/dev/null | ./sito2 | wc -l
3
student@ANT:~/nauty/nauty27r3$ ./geng -c 6 2>/dev/null | ./sito2 | wc -l
6
student@ANT:~/nauty/nauty27r3$ ./geng -c 7 2>/dev/null | ./sito2 | wc -l
7
student@ANT:~/nauty/nauty27r3$ ./geng -c 8 2>/dev/null | ./sito2 | wc -l
22
student@ANT:~/nauty/nauty27r3$ ./geng -c 9 2>/dev/null | ./sito2 | wc -l
24
 
student@ANT:~/nauty/nauty27r3$ time ./geng -c 9 2>/dev/null | ./sito2 | wc -l
24
 
real	0m18,066s
user	0m18,206s
sys	0m0,011s
 
zwierzak@PUTS319:~/nauty2_8_9$ time ./geng -c 9 2>/dev/null | ./sito3 | wc -l
24
 
real    0m11.228s
user    0m11.303s
sys     0m0.019s
 
zwierzak@PUTS319:~/nauty2_8_9$ time ./geng -c 9 2>/dev/null | ./sito5 | wc -l
24
 
real    0m0.745s
user    0m0.809s
sys     0m0.011s
 
zwierzak@PUTS319:~/nauty2_8_9$ time ./geng -c 9 2>/dev/null | ./sito5 | wc -l
24
 
real    0m0.307s
user    0m0.366s
sys     0m0.016s
 
------------------------
 
Wyniki można porównać z 
https://oeis.org/A064731
 
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#define NMAX 20
#define BUFSIZE 1024
 
/**  Funkcja eigensymmatrix jest wzorowana na kodzie w języku Pascal A.Marciniaka */
/** 
@brief Testowanie czy graf jest całkowity 
@param BUFOR - łańcuch tekstowy kodujący graf w formacie graph6 
@return 1 lub 0 
*/
int eigensymmatrix(char * BUFOR)
{
  int i,j,k,k3,k4,L,L1,z;
  long double /* lambda, */ eps,g,h,ma,mn,norm,s,t,u,w;
  int cond;
  long double d[NMAX+1], e[NMAX+1], e2[NMAX+1], Lb[NMAX+1];
  long double x[NMAX+1];
  long double a[NMAX*(NMAX-1)/2 + NMAX + 1];
  int n;
  int bit, poz, poz2;
  bit = 32;
  poz = 1;
  poz2 = 1;
  n = BUFOR[0] - 63;
  a[0] = 0.0;
  for (i = 0; i < n; i++)
   for (j = 0; j<=i; j++)
    {
	if (i==j) {a[poz2++] = 0; }
	else {
         if (bit == 0) { bit = 32;  poz++; }
         if ((BUFOR[poz] - 63) & bit)
                { a[poz2++] = 1; }
               else
                { a[poz2++] = 0; }
         bit = bit >> 1;
	}
    }
 
  int k1 = 1;
  int k2 = n;
  // if ((1<=k1) && (k1<=k2) && (k2<=n))
   {
    i = 0;
    for (L=1;L<=n;L++) { i += L; d[L] = a[i]; /* printf("%Lf ",a[i]); */ }
 
    for (L=n;L>=2;L--)
     {
      i--; j = i; h = a[j]; s = 0;
      for (k=L-2;k>=1;k--) { i--; g = a[i]; s += g*g; }
      i--;
      if (s == 0) { e[L] = h; e2[L] = h*h; a[j] = 0.0; }
       else
        {
          s += h*h; e2[L] = s; g = sqrt(s); if (h>=0.0) g=-g;
          e[L] = g;
          s = 1.0 / (s-h*g);
          a[j] = h - g; h = 0.0; L1 = L - 1; k3 = 1;
          for (j=1;j<=L1;j++)
           {
             k4 = k3; g = 0;
             for (k=1;k<=L1;k++) { g +=a[k4]*a[i+k]; 
                                   if (k<j)  z = 1; else z = k;
                                   k4 += z; }
             k3 += j; g *= s; e[j] = g; h += a[i+j]*g;
           }
          h *= 0.5*s; k3 = 1;
          for (j=1;j<=L1;j++)
           {
             s = a[i+j]; g = e[j]-h*s; e[j] = g;
             for (k=1;k<=j;k++) { a[k3] += -s*e[k]-a[i+k]*g; k3++; }
           }
        }
      h = d[L]; d[L] = a[i+L]; a[i+L] = h;
     }
    h = d[1]; d[1] = a[1]; a[1] = h; e[1] = 0.0; e2[1] = 0.0; s = d[n];
    t = fabs(e[n]); mn = s - t; ma = s + t;
    for (i=n-1;i>=1;i--)
     {
      u = fabs(e[i]); h = t + u; t = u; s = d[i]; u = s - h;
      if (u < mn) mn = u;
      u = s + h;
      if (u > ma) ma = u;
     }
    for (i=1;i<=n;i++) { Lb[i] = mn; x[i] = ma; }
    norm = fabs(mn); s = fabs(ma);
    if (s>norm) norm = s;
    w = ma; /* lambda = norm; */ eps = 7.28e-17*norm;
    for (k=k2;k>=k1;k--)
     {
      /* eps = 7.28e-17*norm; */ s = mn; i = k;
      do {cond = 0; g = Lb[i];
         if (s < g) s = g; else { i--; if (i>=k1) cond = 1; }
      } while (cond);
      g = x[k];
      if (w>g) w = g;
      while (w-s>2.91e-16*(fabs(s)+fabs(w))+eps)
       {
         if (floor(w+10e-5)<s-10e-5) return 0;  // przedział nie zawiera liczby całkowitej
         L1 = 0; g = 1.0; t = 0.5*(s+w);
         for (i=1;i<=n;i++)
          {
            if (g!=0)  g = e2[i] / g; else g = fabs(6.87e15*e[i]);
            g = d[i]-t-g;
            if (g<0) L1++;
          }
         if (L1<k1) { s = t; Lb[k1] = s; }
          else
           { if (L1<k)
               {
                 s = t; Lb[L1+1] = s;
                 if (x[L1]>t) x[L1] = t;
               }
              else w = t;
           }
      } // while
      u = 0.5*(s+w); x[k] = u;
	  if  (!(( ceil(u) - u  < 10e-5 ) || ( u - floor(u) < 10e-5 ))) { return 0; };
    }
  } 
  return 1;
}
 
/** 
@brief Main function: 
* <a href="https://en.cppreference.com/w/cpp/language/main_function">main()</a> 
@param argc
@param argv
@return <a href="https://en.cppreference.com/w/c/program/EXIT_status">EXIT_SUCCESS or EXIT_FAILURE</a>
*/
int main(int argc, char *argv[])
{
 
#ifdef _OPENMP
  //double start; 
  //double end;
  int t = 1;
  char ** bufory = NULL;
 
  if (argc > 1)  t = atoi(argv[1]);
  //start = omp_get_wtime(); 
   bufory = (char **)malloc(t*sizeof(char *));
      //printf("%d\n",t);
   for (int i=0;i<t;i++) bufory[i]=(char*)malloc(BUFSIZE);
 
  char * res = NULL;
  int graphs = 0;
  int i;
  int gid;
 
 do {  
    for (graphs=0;graphs<t;graphs++) {if (fgets(bufory[graphs],BUFSIZE,stdin)==NULL) break;}
    #pragma omp parallel default(none) shared(bufory,graphs) private(gid)
    {
       //#pragma omp for schedule(dynamic,1)
       //#pragma omp for schedule(static,1000)
       #pragma omp for 
	for (gid = 0;gid < graphs;gid++)
	   if (eigensymmatrix(bufory[gid])) {
             // #pragma omp critical 
   	     printf("%s",bufory[gid]);
        }		
    }
 } while(graphs==t);	
 
 for (int i=0;i<t;i++) free(bufory[i]);
 free(bufory);
 //end = omp_get_wtime(); 
  //fprintf(stderr,"| Czas(OpenMP): | %f[sec] |\n",end - start);
#else
  char BUFOR[BUFSIZE];
  // printf("KMAX:%d\n",NMAX*(NMAX-1)/2);
  //clock_t start,stop;
  //start = clock(); 
  while (fgets(BUFOR,BUFSIZE,stdin)) { 
   if (eigensymmatrix(BUFOR)) printf("%s",BUFOR);         
  } // while
  //stop  = clock(); 
  //fprintf(stderr,"| Czas: | %f[sec] |\n",(float)(stop-start)/CLOCKS_PER_SEC);
#endif 
  return EXIT_SUCCESS;
}