#include <math.h>
#include <cmath>
//#include <RcppEigen.h>
#include <RcppArmadillo.h>
//#include <iostream>
//#include <vector>
//#include <Eigen/Core>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::Triplet<float> T; //tripleta para rellenar matrices sparse
typedef Eigen::MappedSparseMatrix<float> MSpMat; //alias para matriz sparse
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd; //alias para vector sparse
using namespace Rcpp;
//using namespace RcppArmadillo;

// [[Rcpp::export]]
float KernelGauss (arma::cube x, arma::cube y, int dimension, double sigJ, double sigd, float distancia)
{
  //Calculo individual del kernel entre dos observaciones
  // input: x observaciones
  //        y observaciones
  // dimension: numero de canales
  //sigJ: sd de la imagen
  //sigd: sd de la posicion
  //distancia: distancia L_2 entre los dos pixeles
  float kernel_value = 0;
  float siJ = (float)sigJ;
  float sid = (float)sigd;
  for (int i = 0; i < dimension; i++ )
  {
    kernel_value +=(float)pow( (float)(x[i] - y[i]),2);
  }
  kernel_value = exp(-(kernel_value/(2*pow(siJ,2)) ) - distancia/(2*pow(sid, 2)));
  return( kernel_value );
}


  
// [[Rcpp::export]]
Eigen::SparseMatrix<float> Kernel_RGB ( arma::cube  M, int h, int w, double distancia, double sigJ, double sigd, int dimension)
{
  // Funcion para calcular  la matriz de similaridades de todos los puntos de una imagen 
  //  input: w: altura en pixeles de la imagen//SINO NO ME SALIAN LAS MULTIPLICACIONES
  //         h: ancho en pixeles de la imagen//SINO NO ME SALIAN LAS MULTIPLICACIONES
  //         distancia: limite de la vecindad del pixel a considerar
  //         M: imagen de 'dimension' dimensiones (numero de canales)
  //         sigJ: sigma para kernel gaussiano, de la imagen
  //         sigd: sigma para kernell gaussiano, de las distancias
  //         dimension: numero de canales
  //output: W, similaridad calculada 
  int i, j, k, l = 0 ;
  float vecindad = distancia;
  float d2 = 0.;
  float calculo = 0.;
  std::vector<T> tripletList; //inicializamos tripleta
  tripletList.reserve((int)floor(h*h*w*w*.1));
  
  //comienza calculo de matriz de similaridades
  for(i=0; i < h; i++)
  {
    for(j=0; j< w; j++ )
    {
      for(k=0; k<h; k++)
      {
        for(l=0; l<w; l++)
        {
          d2 = (float)(pow(i-k,2)+pow(j-l,2));
          if(d2 <= vecindad)
          {
            calculo = KernelGauss( M.tube(j,i), M.tube(l, k), dimension, sigJ, sigd,  vecindad);
            tripletList.push_back(T((j)*h+i, (l)*h+k, calculo ));
          }
        }
      }
    }
  }
  //termina calculo de matriz de similaridades
  Eigen::SparseMatrix<float>  W(h*w, h*w); //construccion de matriz de similaridades
  W.setFromTriplets(tripletList.begin(), tripletList.end());
  return(W);
}



// [[Rcpp::export]]
Eigen::SparseMatrix<float> Kernel_float ( NumericMatrix M, int h, int w, double distancia, double sigJ, double sigd)
{
  // Funcion para calcular  la matriz de similaridades de todos los puntos de una imagen 
  //  input: w: altura en pixeles de la imagen//SINO NO ME SALIAN LAS MULTIPLICACIONES
  //         h: ancho en pixeles de la imagen//SINO NO ME SALIAN LAS MULTIPLICACIONES
  //         distancia: limite de la vecindad del pixel a considerar
  //         M: imagen de 2 dimensiones (EN ESCALA DE GRISES)
  //         sigJ: sigma para kernel gaussiano, de la imagen
  //         sigd: sigma para kernell gaussiano, de las distancias
  //output: W, similaridad calculada, sparse 
  int i, j, k, l = 0 ;
  float vecindad = (float)distancia;
  float d2 = 0.;
  float calculo = 0.;
  std::vector<T> tripletList; //inicializamos tripleta
  tripletList.reserve((int)(h*h*w*w*.1));
  //comienza calculo de matriz de similaridades
  for(i=0; i < h; i++)
  {
    for(j=0; j< w; j++ )
    {
      for(k=0; k<h; k++)
      {
        for(l=0; l<w; l++)
        {
          d2 = (float)(pow(i-k,2)+pow(j-l,2));
          if(d2 <= vecindad)
          {
            calculo = (float) (exp(-pow(M(j,i)-M(l,k),2)/(2*pow(sigJ,2))-d2/(2*pow(sigd,2))));
            tripletList.push_back(T((j)*h+i, (l)*h+k, calculo ));
          }
        }
      }
    }
  }
  //termina calculo de matriz de similaridades
  Eigen::SparseMatrix<float>  W(h*w, h*w); //construccion de matriz de similaridades
  W.setFromTriplets(tripletList.begin(), tripletList.end());
  return(W);
}


