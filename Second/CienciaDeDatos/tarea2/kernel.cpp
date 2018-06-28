#include <Rcpp.h>
//#include <RcppArmadillo.h>
#include <math.h>
#include <cmath>
# define M 1000
using namespace Rcpp;

// [[Rcpp::export]]
double KernelGauss (NumericVector x, NumericVector y, int dimension)
{
    //Calculo individual del kernel entre dos observaciones
        // input: x observaciones
        //        y observaciones
        //double alpha = 0.0013984457316839;
        double alpha = 1.0;
        double kernel_value = 0;
        for (int i = 0; i < dimension; i++ )
        {
          kernel_value += pow( (double)(x[i] - y[i]),2);
        }
        kernel_value = exp(-kernel_value/(2*pow(alpha, 2)));
        
  return( kernel_value );
}

/*R
KernelGauss <- function(x,y, 
                         # alpha = (.05/2)**.5) #kernel gaussiano
                          alpha = 0.0013984457316839 ) #para el caso de vinos
{
  return(exp(-sum((x-y)**2)/(2*alpha**2)))
}
*/

// [[Rcpp::export]]
   NumericMatrix CalculaKernel (NumericMatrix MK, int n , NumericMatrix data, int p)
   {
     // Funcion para calcular  la matriz de kernel de todos los puntos
             //  input: MK matriz inicializada a -1 de tamanio n*n 
             //         n- numero de registros en el dataset
             //output: MK, kernel calculado 
     int i,j = 0  ;
     double calculo=0;
     for(  i = 0; i < n; i++ )
     {       
       for( j= i; j < n; j++) //como el kernel es simetrico solo calculamos la mitad de entradas de la matriz
       {
         calculo = KernelGauss( data.row(i), data.row(j), p ) ;
         MK(i, j) = calculo;
         MK(j, i) =calculo; 
       }
     }
     return(MK);
   }