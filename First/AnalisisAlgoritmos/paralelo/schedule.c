#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>

#ifdef _OPENMP
  #include <omp.h>
  #define TRUE  1
  #define FALSE 0
#else
  #define omp_get_thread_num() 0
#endif

int main()
{
   int i, j, n = 9;

#ifdef _OPENMP
   (void) omp_set_dynamic(FALSE);
   if (omp_get_dynamic()) {printf("Advertencia: se ha hecho el ajuste dinamico de hilos\n");}
   (void) omp_set_num_threads(4);
#endif

#pragma omp parallel for default(none) schedule(runtime) \
        private(i,j) shared(n)
   for (i=0; i<n; i++)
   {


      for (j=0; j<i; j++)
          Sleep(2000);
          printf("Iteracion %d ejecutada por el hilo %d\n",
          i, omp_get_thread_num());
   } // Final del for paralelo

   return(0);
}
