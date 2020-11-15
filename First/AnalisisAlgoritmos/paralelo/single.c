#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
  #include <omp.h>
  #define TRUE  1
  #define FALSE 0
#else
  #define omp_get_thread_num() 0
#endif

int main()
{

#ifdef _OPENMP
   (void) omp_set_dynamic(FALSE);
   if (omp_get_dynamic()) {printf("Warning: dynamic adjustment of threads has been set\n");}
   (void) omp_set_num_threads(4);
#endif

   int n = 9;
   int i, a, b[9] = {-1};

#pragma omp parallel shared(a,b, i)
{
   #pragma omp single
   {
      a = omp_get_thread_num();
      printf("Constructo single ejecutado por el hilo %d\n",
             omp_get_thread_num());
   }
   #pragma omp for
   for (i=0; i<n; i++)
       b[i++] += omp_get_thread_num();

} // Final de la region paralela

   printf("Despues de la region paralela:\n");
   for (i=0; i<n; i++)
       printf("b[%d] = %d\n",i,b[i]);

   return(0);
}
