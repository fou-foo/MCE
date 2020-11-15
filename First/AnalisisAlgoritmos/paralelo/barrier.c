#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <strings.h>
#include <Windows.h>

#ifdef _OPENMP
  #include <omp.h>
  #define TRUE  1
  #define FALSE 0
#else
  #define omp_get_thread_num() 0
  #define omp_get_num_threads() 1
#endif

void print_time(int TID, char *comment);

int main()
{
   int TID;
   int i, n = 10;
   int a[n], b[n], ref[n];

#ifdef _OPENMP
   (void) omp_set_dynamic(FALSE);
   if (omp_get_dynamic()) {printf("Advertencia: se ha hecho el ajuste dinamico de hilos\n");}
   (void) omp_set_num_threads(4);
#endif

   #pragma omp parallel private(TID)
   {

      TID = omp_get_thread_num();

      if ( TID < omp_get_num_threads()/2 ) Sleep(1000);

      (void) print_time(TID,"antes");

      #pragma omp barrier

      (void) print_time(TID,"despues ");

   } //Final de la region paralela

   for (i=0; i<n; i++)
   {
       b[i] = 2*(i+1);
       ref[i] = i + b[i];
   }

   #pragma omp parallel private(i) shared(n,a,b)
   {
      #pragma omp for schedule(dynamic,1) nowait
      for (i=0; i<n; i++)
          a[i] = i;

      #pragma omp barrier

      #pragma omp for schedule(dynamic,1) nowait
      for (i=0; i<n; i++)
          a[i] += b[i];
   } //Final de la region paralela

   printf("Despues de la region paralela\n");
   for (i=0; i<n; i++)
       printf(" a[%3d] = %6d ref[%3d] = %6d\n",i,a[i],i,ref[i]);

   return(0);
}

void print_time(int TID, char *comment)
{
   time_t tp;
   char buffer[26], mytime[9];

   (void) time(&tp);
   strcpy(&buffer[0],ctime(&tp));

   strncpy(&mytime[0],&buffer[11],8);
   mytime[8]='\0';

   printf("Hilo %d %s de la barrera en %s\n",TID,comment,&mytime[0]);

   return;
}
