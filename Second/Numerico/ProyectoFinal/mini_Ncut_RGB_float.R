setwd('/home/fou/Desktop/MCE/Second/Numerico/ProyectoFinal')
library(Rcpp) #libreria para codigo c++
library(RSpectra) #libreria para lanczos
library(RcppEigen) #libreria para codigo c++
library(imager) #libreria para leer imagenes
library(Matrix)
sourceCpp('W_RGB_float.cpp') #compilamos el programa en C++
t1 <- Sys.time() #medimos tiempo de ejecucion
dir()
imagen <- load.image('001.jpg')
dim(imagen)
plot(imagen) #visualizamos la imagen
#############preprosesamiento  interpolacion 
imagen2 <- imagen
imagen2 <- resize(im = imagen, size_x = 3001, size_y = 2257, size_z = 1, size_c = 3 )
#imagen2 <- resize_halfXY(imagen)
#imagen2 <- resize_halfXY(imagen2)
#imagen2 <- resize_halfXY(imagen2)

# Estandarizacion RGB canal por canal: aumentar contraste
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
## definimos parametros
(h <- 94)
(w <- 120 )
porcentaje <- .1
edges <- ((h*w)*(h*w-1)/2)*(porcentaje) #segun el paper se pueden remover hasta el 90% de las aristas deberia de ser ((h*w)*(h*w+1)/2)*.1
edges.sugerido <- edges/(h*w) #promedio de aristas por nodo
cuadrado <- edges.sugerido**.5 # vamos a fijar esta cantidad
cuadrado <- round(cuadrado) +1 
sigJ <- 0.05 #ver paper
sigd <- 10#ver paper
r2 <- cuadrado**2
eje.x <- 0:24*120+1
eje.y <- 0:23*94+1
copia <- gray.imagen[,,1,1:3]
class(copia)
dim(copia)
for( paso.x in eje.x)
{
    for(paso.y in eje.y)
    {
        gc()
        t1 <- Sys.time()
        M1 <-  copia[paso.x:(paso.x+119), paso.y:(paso.y+93), 1:3 ]
        dim(M1)
        plot(as.cimg(M1))
        #scan()
        W <- Kernel_RGB(M1, h, w, r2, sigJ, sigd, dimension=3)
        W <- as(W, "sparseMatrix")
        print('tick1')
        remove(M1)
        gc()
        #isSymmetric(W)
        d <- Matrix::colSums(W) #obtenemos suma por columnas
        D_medio <- Matrix::.sparseDiagonal(n = h*w, x = d**(-.5) ) #calculamos la matriz D^{-1/2} para el problema de valores propios generalizado
        W <- D_medio%*%(Matrix::.sparseDiagonal(n = h*w, x = d ) -W)%*%D_medio 
        print('tick2')
        Z <- eigs_sym(W, k=2, which='LM', sigma = 0) #usamos implicited restarted lanczos
        Z$values #visualizamos los dos valores propios mas pequenios
        remove(W) #ahorramos RAM
        remove(d)
        gc()
        print('tick3')
        #####################
        Y1 <- D_medio%*%(Z$vectors[,1 ]) # Usar los vectores propios que se encontraron para segmentar
        hist(as.matrix(Y1))
        #Y2 <- D_medio%*%Z$vectors[,2]
        #hist(as.matrix(Y2))
        set.seed(0)
        mascara <- Y1>0
        mascara <- matrix(mascara, ncol = h, byrow = TRUE)
        segmentacion <- mascara
        table(segmentacion)
        imagen.segmentacion <- as.cimg(round(segmentacion,1))
        plot(imagen.segmentacion)
        Aplica.Mascara <- function( x){
            M1 <- copia[paso.x:(paso.x+119), paso.y:(paso.y+93) ,x]
            M1 <- as.matrix(M1)*mascara
            copia[paso.x:(paso.x+119), paso.y:(paso.y+93) ,x] <<- M1
            return((copia[ , ,x]))
        }
        copia <- lapply(1:3, FUN = Aplica.Mascara)
        library(abind)
        copia <- abind(copia , along = 3)
        imagen.final <- as.cimg(copia)
        plot(imagen.final)
        set.seed(0)
        #scan()
        t1 <- Sys.time()-t1
        print(t1)
        print('paso x')
        print(paso.x)
        print('paso y')
        print(paso.y)
        #scan()
        save(copia, file = 'Ncut.RData')
        plot(as.cimg(copia))
        #foo <- scan()
        gc()
    }
}
plot(as.cimg(copia))
