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
   int i, a, n = 5;

#ifdef _OPENMP
   (void) omp_set_dynamic(FALSE);
   if (omp_get_dynamic()) {printf("Advertencia: se ha hecho el ajuste dinamico de hilos\n");}
   (void) omp_set_num_threads(3);
#endif

    #pragma omp parallel shared(i)
   {
       #pragma omp parallel for  lastprivate(a)
   for (i=0; i<n; i++)
   {
       a = i+1;
       printf("El hilo %d tiene un valor de a = %d para i = %d\n",
              omp_get_thread_num(),a,i);
   } // Final del for paralelo
   }
   printf("Valor de a después del for paralelo: a = %d\n",a);

   return(0);
}
