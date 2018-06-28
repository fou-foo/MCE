setwd('C:\\Users\\fou-f\\Desktop\\MCE\\Second\\CienciaDeDatos\\tarea2')
strings <- read.csv('train_stock.csv', stringsAsFactors = FALSE)
## veamos la longitud en palabras de los clasificados como iguales
iguales <- subset(strings, same_security==TRUE) #sleccionamos solo los etiquetados como semejantes
x1.len <- lapply(FUN= function(i){
    n <- strsplit(as.character(iguales[i,'description_x'])," " )
    n <- unlist(n)
    n.len <- length(n)
    n.len
}, X =(1:dim(iguales)[1])) # checamos la longitud en palabras de la descripcion x
x2.len <- lapply(FUN= function(i){
  n <- strsplit(as.character(iguales[i,'description_y'])," " )
  n <- unlist(n)
  n.len <- length(n)
  n.len
}, X =(1:dim(iguales)[1])) # checamos la longitud en palabras de la descripcion y
x1.len <- unlist(x1.len)
x2.len <- unlist(x2.len)
(table(x1.len, x2.len))# 4 es el mayor 
## veamos la longitud en palabras de los clasificados como diferentes
diferentes <- subset(strings, same_security==FALSE)
y1.len <- lapply(FUN= function(i){  # checamos la longitud en palabras de la descripcion x
  n <- strsplit(as.character(diferentes[i,'description_x'])," " )
  n <- unlist(n)
  n.len <- length(n)
  n.len
}, X =1:dim(diferentes)[1])
y2.len <- lapply(FUN= function(i){ # checamos la longitud en palabras de la descripción y
  n <- strsplit(as.character(diferentes[i,'description_y'])," " )
  n <- unlist(n)
  n.len <- length(n)
  n.len
}, X =1:dim(diferentes)[1])
y1.len <- unlist(y1.len)
y2.len <- unlist(y2.len)
table(y1.len, y2.len) #5 es el mayor
library(kernlab)
# consideramos 4 familias de kernels, con mayores esperanzas en el spectrum de longitud fija y lambda=1
#encontramos de mayor utilidad los kernels normalizados
spectrum.kernel <- stringdot(type="spectrum",length=3,
                             normalized=TRUE, lambda=1)
constante.kernel <- stringdot(type='constant', normalized = FALSE)
boundrange.kernel <- stringdot(type ='boundrange', length = 10,
                               normalized = TRUE, lambda = 1)
exponential.kernel <- stringdot(type='exponential', 
                                lambda = .5, normalized = TRUE)

k <- spectrum.kernel
#construimos la primer matriz de similaridades, con todos los datos del archivo
K.Gram <- kernelMatrix(kernel=k ,
                       x=strings$description_x, 
                       y=strings$description_y)
#realizamos Kernel - pca considerando un kernel gaussiano con parámetro igual a la
#unidad, porque la matriz esta de similaridades esta centrada
kernel.pca <- kpca(K.Gram, kernel = "rbfdot", kpar = list(sigma = 1)) 
# revisamos la varianza explicada por las 10 primeras componentes principales
val.pro <- kernel.pca@eig
plot(cumsum(kernel.pca@eig)[1:10]/sum(kernel.pca@eig))
#guardamos los datos rotados: NOTA no tiene porque ser cuadrada la 
# matriz de vectores propios, pues depende del rango de la matriz de 
# similaridades
rotacion <- kernel.pca@rotated
rotacion <- as.data.frame(rotacion)
rotacion$etiqueta <- factor(strings$same_security)
library(ggplot2) #graficamos la proyeccion en las dos primeras componentes
ggplot(rotacion, aes(x =V1, y = V2, color=etiqueta, alpha=.1))+ geom_point()+
  ggtitle('Kernel: spectrum, n=3, lambda=1') + theme_minimal() + xlab('Primer kernel-componente principal')+
  ylab('Segunda kernel-componente principal')+ guides(alpha = FALSE)
###########inciso b 
set.seed(0) #fijamos la semilla para comparar la pseudo muestra aleatoria
K.Gram.train <- kernelMatrix(kernel=k , #construimos la matriz de similaridades 
                       x=strings$description_x) #considerando solo las etiquetas 'description_x'
kernel.pca.train <- kpca(K.Gram.train, kernel = "rbfdot", kpar = list(sigma = 1)) #kernel pca igual que en el caso anterior
index <- sample(1:dim(strings)[1], dim(strings)[1]  ) #indice de la muestra aleatoria
componentes.prin.train <- kernel.pca.train@pcv #vectores propios de kernel-pca
test <- strings[index, ] #conjunto de test
K.Gram.test <- kernelMatrix(kernel=k ,y = strings$description_x,
                             x=test$description_y) #calculamos la similaridad entre el 
                #conjunto de train y el de test
res <- predict(kernel.pca.train, K.Gram.test) #rotamos el conjunto test
    #se procede a calculo explicito de las distancias de cada observacion
    #en el conjunto test contra todas las observaciones rotadas en el conjunto de train
foo <- matrix(rep(-1, dim(res)[1]*dim(strings)[1]), ncol = dim(strings)[1])
for (i in 1:dim(foo)[1])
{
  for(j in 1:dim(foo)[2])
  {
    foo[i, j] <- sum((res[i,]-kernel.pca.train@rotated[j,])**2)**.5
  }
}
mas.proximo <- apply(foo, 1, which.min ) #determinamos la observacion mas cercana del conjunto train
train <- kernel.pca.train@rotated[,1:2]
train <- as.data.frame(train)
test$res <- strings$description_x[mas.proximo] #agregamos al test para comparar
test$resul <- test$description_x == test$res #calificamos prediccion vs etiqueta
round(table(test$same_security,test$resul)/sum(table(test$same_security,test$resul)),2) #resumen 
table(strings$same_security)/dim(strings)[1]
res <- as.data.frame(res)
y_hat <- train[mas.proximo,]
vis <- cbind(y_hat, res)
colores <- c('purple', 'black', 'orange', 'red', 'turquoise', 'blue4',
             'green4', 'salmon')
vis$index <- rep(colores, 3)
ggplot(train, aes(x =V1, y = V2, color=I('lightblue'), alpha=1))+ geom_point()+
  ggtitle('Train') + theme_minimal() + xlab('Primer kernel-componente principal')+
  ylab('Segunda kernel-componente principal')+ guides(alpha = FALSE) +
  geom_point(data=res[, 1:2], aes(x =V1, y = V2, 
                                  color=(vis$index), alpha =1))+
  geom_point(data=y_hat[, 1:2], aes(x =V1, y = V2, 
                                    color=(vis$index), alpha =1))+
  xlim(c(-26, 13))+ylim(c(-20,15))

#spectrum 1 lambda =1: .224 .2009
#spectrum 2 lambda =1: .224 .2289
#spectrum 3 lambda =1: .219 .2336
#spectrum 4 lambda =1: .219 .2289  
#boundrange 1 : .2242 .2009
#boundrange 2 : .219 .214
#boundrange 3 : .219 .219
#boundrange 5 : .2242 .228
#boundrange 10 : .2242 .2336
#constant : .2336 .112
#exponential 10 : .224 .2102 
#exponential 20 : .224 .2009 
#exponential .5 :  .233 .163