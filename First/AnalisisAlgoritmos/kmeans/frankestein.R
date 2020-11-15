library(Rcpp)
library(RcppArmadillo)
sourceCpp("C:\\Users\\fou-f\\Desktop\\MCE\\First\\AnalisisAlgoritmos\\kmeans\\kmeansC.cpp")
source("C:\\Users\\fou-f\\Desktop\\MCE\\First\\AnalisisAlgoritmos\\kmeans\\kmeans_r.R")

a <- kmeans(iris[, 1:4], 3)
iteraciones <- a$iter
a <- CfooKmeans(data_m, iteraciones, distancias,  centroides, n, columnas_vector, columnas, k )

microbenchmark( 
  
    CFoomeans =CfooKmeans(data_m, iteraciones, distancias,  centroides, n, columnas_vector, columnas, k ),
    R = kmeans(data[,1:4], iter.max = 10, centers = k , algorithm = c("MacQueen")),
   Foomeans = kmeans.foo(data, k, iteraciones, columnas.vector= c(1,2,3,4) )  ,
            times = 1000L          )
               
               
               
               

