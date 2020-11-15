
#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
  #include <omp.h>
  #define TRUE  1
  #define FALSE 0
#endif

int main()
{


#ifdef _OPENMP
   (void) omp_set_dynamic(FALSE);
   if (omp_get_dynamic()) {printf("Warning: dynamic adjustment of threads has been set\n");}
   (void) omp_set_num_threads(4);
#endif

   int i, n = 7;
   int a[n];

   for (i=0; i<n; i++)
      a[i] = i+1;
for (i=0; i<n; i++)
      printf("a[%d] = %d\n",i,a[i]);
   #pragma omp parallel for shared(a)
   for (i=0; i<n; i++)
   {
       a[i] += i;
   } // Final del for paralelo

   printf("En el programa principal despues del for paralelo para:\n");
   for (i=0; i<n; i++)
      printf("a[%d] = %d\n",i,a[i]);

   return(0);
}
