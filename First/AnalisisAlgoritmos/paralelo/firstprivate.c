#include <stdio.h>
#include <stdlib.h>

#define TRUE  1
#define FALSE 0

#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
  #define omp_get_num_threads() 1
#endif

int main()
{


#ifdef _OPENMP
   (void) omp_set_dynamic(FALSE);
   if (omp_get_dynamic()) {printf("Advertencia: se ha hecho el ajuste dinamico de hilos\n");}
   (void) omp_set_num_threads(3);
#endif

   int *a;
   int n = 2, nthreads, vlen, indx, offset = 4, i, TID;
   int failed;

   indx = offset;

   // Establecer parametros y asignar memoria

   #pragma omp parallel firstprivate(indx) shared(a,n,nthreads,failed)
   {
     #pragma omp single
     {
        nthreads = omp_get_num_threads();
        vlen = indx + n*nthreads;
        if ( (a = (int *) malloc(vlen*sizeof(int))) == NULL )
           failed = TRUE;
        else
           failed = FALSE;
     }
   } // Final de la region parela

   if ( failed == TRUE ) {
      printf("Error: la asignaciond de memoria fallo cuando vlen = %d\n",vlen);
      return(-1);
   }
   else
   {
      printf("Diagnosticos:\n");
      printf("nthreads = %d\n",nthreads);
      printf("indx     = %d\n",indx);
      printf("n        = %d\n",n);
      printf("vlen     = %d\n",vlen);
   }

   for(i=0; i<vlen; i++)
        a[i] = -i-1;

   // Cada hilo accesa a un arreglo mediante la variable indice indx

   printf("La longitud del segmento por hilo es %d\n",n);
   printf("El offset para el vector a es %d\n",indx);
   #pragma omp parallel default(none) firstprivate(indx) \
           private(i,TID) shared(n,a)
   {
      TID = omp_get_thread_num();
      indx += n*TID;
      for(i=indx; i<indx+n; i++)
         a[i] = TID + 1;
   } // Final de la region paralela

   printf("Despues de la region paralela:\n");
   for (i=0; i<vlen; i++)
      printf("a[%d] = %d\n",i,a[i]);

   free(a);

   return(0);
}
