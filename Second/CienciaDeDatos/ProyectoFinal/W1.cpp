#include <math.h>
#include <cmath>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <iostream>
// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::Triplet<float> T; //tripleta para rellenar matrices sparse
typedef Eigen::MappedSparseMatrix<float> MSpMat; //alias para matriz sparse
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd; //alias para vector sparse
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix Kernel_float ( NumericMatrix M, int h, int w, double distancia, double sigJ, double sigd)
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
//  std::vector<T> tripletList; //inicializamos tripleta
//  tripletList.reserve((int)floor(h*h*w*w*.1));
  //comienza calculo de matriz de similaridades
  for(i=0; i < h; i++)
  {
    printf("i :  %d  \n ", i);
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
  //          tripletList.push_back(T((j)*h+i, (l)*h+k, calculo ));
          }
        }
      }
    }
  }
  //termina calculo de matriz de similaridades
//  Eigen::SparseMatrix<float>  W(h*w, h*w); //construccion de matriz de similaridades
//  W.setFromTriplets(tripletList.begin(), tripletList.end());
  return(M);
}