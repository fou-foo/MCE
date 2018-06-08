setwd('C:\\Users\\fou-f\\Desktop\\MCE\\Second\\AnalisisNumericoYOptimizacion\\Miniproyecto')
library(Rcpp) #libreria para codigo c++
library(RcppEigen) #libreria para codigo c++
library(RSpectra) #libreria para lanczos
library(imager) #libreria para leer imagenes
library(Matrix)
sourceCpp('W_RGB_float.cpp') #compilamos el programa en C++
t1 <- Sys.time() #medimos tiempo de ejecucion
dir()
imagen <- load.image('marco.jpg')
dim(imagen)
plot(imagen) #visualizamos la imagen
prod(dim(imagen)[1:2]/12.1) -160*120
#imagen2 <- imagen
#############preprosesamiento  interpolacion 
imagen2 <- imagen
#imagen2 <- resize(im = imagen, size_x = 166, size_y = 110, size_z = 1, size_c = 3 )
#imagen2 <- resize_halfXY(imagen)
#imagen2 <- resize_halfXY(imagen2)
#imagen2 <- resize_halfXY(imagen2)

#standarizacion RGB canal por canal
imagen2 <- imagen2[, , 1, 1:3] 
estandariza <- function(canal){
  M1 <- imagen2[,,canal]
  M1 <- (M1-min(M1))/(max(M1)-min(M1))
return((M1))
}
M <- lapply(1:3, FUN = estandariza)
library(abind)
M <- abind(M[[1]], M[[2]],M[[3]], along = 3)
gray.imagen <- as.cimg(M)
plot(gray.imagen) #imagen normalizada
##
sig <- 1
siz <- 2*round(3*sig) + 1
#decidimos no usar filtros pasabajas
(h <- dim(M)[2])
(w <- dim(M)[1] )
## definimos parametros
porcentaje <- .0001
edges <- ((h*w)*(h*w-1))*(porcentaje) #segun el paper se pueden remover hasta el 90% de las aristas deberia de ser ((h*w)*(h*w+1)/2)*.1
edges.sugerido <- edges/(h*w) #promedio de aristas por nodo
cuadrado <- edges.sugerido**.5 # vamos a fijar esta cantidad
cuadrado <- round(cuadrado) +1 
sigJ <- 0.05 #ver paper
sigd <- 10#ver paper
#r2 <- 25.*25. #se pueden quitar hasta el 90% de las aristas segun el paper
r2 <- cuadrado**2
dim(M)
W <- Kernel_RGB(M, h, w, r2, sigJ, sigd, dimension=3)
W <- as(W, "sparseMatrix")
print('tick1')
remove(M)
gc()
#isSymmetric(W)
d <- Matrix::colSums(W) #obtenemos suma por columnas
D_medio <- Matrix::.sparseDiagonal(n = h*w, x = d**(-.5) ) #calculamos la matriz D^{-1/2} para el problema de valores propios generalizado
W <- D_medio%*%(Matrix::.sparseDiagonal(n = h*w, x = d ) -W)%*%D_medio 
print('tick2')
Z <- eigs_sym(W, k=3, which='LM', sigma = 0) #usamos lanczos
Z$values #visualizamos los tres valores propios mas pequenios
remove(W) #ahorramos RAM
remove(d)
gc()
print('tick3')
#####################
Y1 <- D_medio%*%(Z$vectors[,1]) #MI GRAN DUDA ERA ESTA, como usar los vectores propios que se encontraron para segmentar
hist(as.matrix(Y1))
Y2 <- D_medio%*%Z$vectors[,2]
hist(as.matrix(Y2))
set.seed(0)
kmeans <- kmeans(unlist(Y2), centers = 2, nstart = 50)
print(table(kmeans$cluster))
mascara <- matrix(kmeans$cluster-1, ncol = h, byrow = TRUE)
segmentacion <- mascara
table(segmentacion)
imagen.segmentacion <- as.cimg(round(segmentacion,1))
Aplica.Mascara <- function( canal){
  M1 <- imagen2[,,canal]
  M1 <- as.matrix(M1)*mascara
  return((M1))
}
M <- lapply( FUN = Aplica.Mascara, X=1:3)
library(abind)
M <- abind(M , along = 3)
imagen.final <- as.cimg(M)
plot(imagen.final)#CON MADRE
set.seed(0)
data.a.segmentar <- cbind(Y1,Y2 )
kmeans2 <- kmeans(data.a.segmentar, centers = 3, nstart = 50)
kmeans2$cluster <- kmeans2$cluster-1
kmeans2$cluster <- kmeans2$cluster/2
print(table(kmeans2$cluster))
mascara <- matrix(kmeans2$cluster, ncol = h, byrow = TRUE)
segmentacion <- mascara
M <- lapply( FUN = Aplica.Mascara, X=1:3)
library(abind)
M <- abind(M , along = 3)
imagen.final <- as.cimg(M)
plot(imagen.final)#CON MADRE
t1 <- Sys.time() -t1
print(t1)
gc()

