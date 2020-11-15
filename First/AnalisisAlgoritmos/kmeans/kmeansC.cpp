#include <Rcpp.h>
//#include <RcppArmadillo.h>
#include <math.h>
#include <cmath>
# define M 1000
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector rngCppScalar() {
  NumericVector x(4);
  x[0] = 1;
  x[1] = 2;
  x[2] = 3;
  x[3] = 4;
  return(x);
}



/*** R
k <- 3    #Numero de kluster deseados  test:3, 4
set.seed(0)
columnas_vector <-  rngCppScalar()
data <- iris                          #Datos con lo que se trabaja
n <- dim(data)[1]                     #Numero de registros en el dataset
set.seed(1)
muestra <- sample(1:n, k)
  */

// [[Rcpp::export]]
NumericMatrix  Ccentroides(int ii, int jj) {
  NumericMatrix x(ii, jj);    //inicializacion de los centroides a cero
  int i,j = 0;
  for (i=0; i<ii; i++)
    for (j=0; j<jj; j++)
      x(i,j) = 0;
  return(x);
}

/*** R
columnas <- 4
centroides <-  Ccentroides(k,columnas)
distancias <-Ccentroides(n,columnas)
*/

// [[Rcpp::export]]
NumericMatrix Cinicializacion (int k, int n, NumericVector  columnas_vector, DataFrame datos, NumericMatrix centroides,NumericVector muestra )
{       //  input: k-numero de centroides int
        //         n- numero de registros en el dataset
        //         columnas-numero de columnas con las que se trabaja,
        //                   y sus indices columnas.vector
        //         data set original
        //  output: >inicializa los centroides a valores en la muestra
  int i,j;
  for (i = 0; i< k ; i++)
  { 
    //NumericVector y(4);
    for ( j=0; j<4;j++) { //checar indices
      NumericVector column = datos[j] ;
      centroides(i, j ) = column[i] ;
    }
  }
  return(centroides);
}




/*** R
centroides <-  Cinicializacion ( k, n, columnas_vector, data, centroides, muestra )
  data_m <- as.matrix(data[,1:4])
  */

// [[Rcpp::export]]
double Cl2 (NumericVector centroide, NumericVector observacion, int columnas)
{
    //Calculo individual de la distancia entre dos vectores
        // input: centroide-uno de los centroides
        //        observacion-uno de los registros del dataset
        double suma = 0.;
        for (int i = 0; i < columnas; i++ )
        {
          suma += pow( (double)(centroide[i]-observacion[i]),2);
        }
        suma = sqrt( (double)suma);
        
  return( suma );
}

 // [[Rcpp::export]]
   NumericMatrix CM_l2 (NumericMatrix distancias, NumericMatrix centroides, int n , int k, NumericMatrix data_m, NumericVector columnas_vector, int columnas)
   {
     // Funcion para actualizar las distancias a los puntos
             //  input: k-numero de centroides int
             //         n- numero de registros en el dataset
             //  input: centroides
             //         data set original
             //         indices columnas.vector
             //output: distancias actualizadas 
             int i,j = 0  ;
     for(  i = 0; i<n; i++ )
     {
       
       for( j=0; j < k; j++)
       {
         distancias(i, j) = Cl2(data_m.row(i), centroides.row(j), columnas) ; 
       }
     }
     return(distancias);
   }

// [[Rcpp::export]]
NumericMatrix Cetiquetamiento (NumericMatrix distancias, int n, int k )
{
    // input : n - numero de observaciones
    //      : k - numero de observaciones
    //# distancias : distancias actualizadas
  
  for( int i = 0; i<n ; i++ )
  {
    distancias(i,k)= M;
    double mini = min( distancias.row(i ) ) ;
    for(int w = 0; w < 4; w++)
    {
      if(distancias(i, w) == mini)
      {
        distancias(i, 3) = w+1;
      }
    }
  }
  return(distancias);
}



/*** R
actualizacion.centroides <- function(data, distancias, k, centroides, columnas.vector )
{
  for(i in 1:k)
  {
    centroides[i,] <- media(data, distancias, i, k, columnas.vector)
  }
  return(centroides)
}


media <- function(data, distancias, centro, k, columnas.vector)
{
# Actualizacion de centroide en cada cluster
  
  subset <- data[ distancias[, k+1] == centro,   columnas.vector ]
  suma <- rep(0, length(columnas.vector))
  for(ii in 1:dim(subset)[1])
  {
    for(jj in 1:length(columnas.vector))
      suma[jj ]<- suma[jj]+ subset[ii,jj] 
  }
  return(suma/dim(subset)[1])
}
*/


// [[Rcpp::export]]
NumericMatrix CfooKmeans(NumericMatrix data_m, int iteraciones, NumericMatrix distancias, NumericMatrix centroides, int n, NumericVector columnas_vector, int columnas, int k )
{
  for( int i =1; i <=iteraciones; i++)
  {
    distancias =  CM_l2( distancias,  centroides,  n , k,  data_m, columnas_vector ,  columnas);
    //particion de voronoi
    distancias = Cetiquetamiento(distancias, n, k ) ;        //asignacion de kluster;
    /**R
    centroides = actualizacion.centroides ( data_m, distancias, k, centroides, columnas_vector );
     */
  }
  return(distancias );
}
