#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
  #include <omp.h>
  #define TRUE  1
  #define FALSE 0
#else
  #define omp_get_thread_num() 0
  #define omp_get_num_threads() 1
#endif

int main()
{

#ifdef _OPENMP
   (void) omp_set_dynamic(FALSE);
   if (omp_get_dynamic()) {printf("Warning: dynamic adjustment of threads has been set\n");}
   (void) omp_set_num_threads(4);
#endif

   int i, n = 9;
   int a[9], b[9];

#pragma omp parallel default(none) shared(n,a,b, i)
{
   #pragma omp single
      printf("Primer ciclo for: el numero de hilos es %d\n",
             omp_get_num_threads());

   #pragma omp for
   for (i=0; i<n; i++)
   {
      printf("El hilo %d ejecuta en el ciclo la iteracion %d\n",
             omp_get_thread_num(),i);
      a[i] = i;
   }

   #pragma omp single
      printf("Segundo ciclo for: el numero de hilos %d\n",
             omp_get_num_threads());

   #pragma omp for
   for (i=0; i<n; i++)
   {
      printf("El hilo  %d ejecuta en el ciclo la iteracion %d\n",
             omp_get_thread_num(),i);
      b[i] = 2 * a[i];
   }
} // Final de la region paralela

   return(0);
}
