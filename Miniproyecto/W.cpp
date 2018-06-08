#include <math.h>
#include <cmath>
//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
//using namespace RcppArmadillo;

// [[Rcpp::export]]
double KernelGauss (arma::cube x, arma::cube y, int dimension, double sigJ, double sigd, double distancia)
{
  //Calculo individual del kernel entre dos observaciones
  // input: x observaciones
  //        y observaciones
  // dimension: numero de canales
  //sigJ: sd de la imagen
  //sigd: sd de la posicion
  //distancia: distancia L_2 entre los dos pixeles
  double kernel_value = 0;
  for (int i = 0; i < dimension; i++ )
  {
    kernel_value += pow( (double)(x[i] - y[i]),2);
  }
  kernel_value = exp(-(kernel_value/(2*pow(sigJ,2)) ) - distancia/(2*pow(sigd, 2)));
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
NumericMatrix CalculaKernel (NumericMatrix W , arma::cube  M, int h, int w, double distancia, double sigJ, double sigd, int dimension)
{
  // Funcion para calcular  la matriz de similaridades de todos los puntos de una imagen 
  //  input: W: matriz inicializada a 0 de tamanio (h*w)^2 
  //         w: altura en pixeles de la imagen//SINO NO ME SALIAN LAS MULTIPLICACIONES
  //         h: ancho en pixeles de la imagen//SINO NO ME SALIAN LAS MULTIPLICACIONES
  //         distancia: limite de la vecindad del pixel a considerar
  //         M: imagen de 'dimension' dimensiones (numero de canales)
  //         sigJ: sigma para kernel gaussiano, de la imagen
  //         sigd: sigma para kernell gaussiano, de las distancias
  //         dimension: numero de canales
  //output: W, similaridad calculada 
  int i, j, k, l = 0 ;
  double vecindad = distancia;
  double d2 = 0.;
  //arma::mat x(1, dimension, arma::fill::zeros);
  //arma::mat y(1, 3, arma::fill::zeros);
  
  for(i=0; i < h; i++)
  {
    for(j=0; j< w; j++ )
    {
      for(k=0; k<h; k++)
      {
        for(l=0; l<w; l++)
        {
          d2 = pow(i-k,2)+pow(j-l,2);
          if(d2 <= vecindad)
          {
            W( (j)*h+i,(l)*h+k)= KernelGauss( M.tube(j,i), M.tube(l, k), dimension, sigJ, sigd,  distancia);

          }
        }
      }
    }
  }
  return(W);
}





// [[Rcpp::export]]
NumericMatrix Kernel (NumericMatrix W , NumericMatrix M, int h, int w, double distancia, double sigJ, double sigd)
{
  // Funcion para calcular  la matriz de similaridades de todos los puntos de una imagen 
  //  input: W: matriz inicializada a 0 de tamanio (h*w)^2 
  //         w: altura en pixeles de la imagen//SINO NO ME SALIAN LAS MULTIPLICACIONES
  //         h: ancho en pixeles de la imagen//SINO NO ME SALIAN LAS MULTIPLICACIONES
  //         distancia: limite de la vecindad del pixel a considerar
  //         M: imagen de 2 dimensiones (EN ESCALA DE GRISES)
  //         sigJ: sigma para kernel gaussiano, de la imagen
  //         sigd: sigma para kernell gaussiano, de las distancias
  //output: W, similaridad calculada 
  int i, j, k, l = 0 ;
  double vecindad = distancia;
  double d2 = 0.;
  for(i=0; i < h; i++)
  {
    for(j=0; j< w; j++ )
    {
      for(k=0; k<h; k++)
      {
        for(l=0; l<w; l++)
        {
          d2 = pow(i-k,2)+pow(j-l,2);
          if(d2 <= vecindad)
          {
            W( (j)*h+i, (l)*h+k) = exp(-pow(M(j,i)-M(l,k),2)/(2*pow(sigJ,2))-d2/(2*pow(sigd,2)));
          }
        }
      }
    }
  }
  return(W);
}