setwd('C:\\Users\\fou-f\\Desktop\\MCE\\Second\\CienciaDeDatos\\tarea2\\data_fruits_tarea')
imagenes <- dir()
library(imager)
L2 <- function(x,y) #funcion para medir distancia entre pixeles
{
  # x (vector-numeric): pixel de la imagen
  # y (vector-numeric): pixel con valores (1,1,1)
  return(sqrt(sum((x-y)**2)))
}
index <- length(imagenes)
medianas <- matrix(rep(-1,3*index), nrow = index) #se reserva espacio para guardar las medianas de los incisos 1 a 3
##########inicia etapa de preprosemiento la cual consiste en trabajar 
##########solo con pixeles que difieren en una pequena distancia del pixel blanco
for(i in 1:index)
{
  fruta <- load.image(file=imagenes[i])
  print(imagenes[i])
  copia <- fruta
  # for(j in 1:100)#para cada pixel se mide su distancia con el pixel blanco
  # {
  #   for(k in 1:100)
  #   {
  #     if(  L2(copia[j,k,,1:3], rep(1.,3)) <=0.02037707)
  #     {
  #       #aproximadamente si el pixel dista en menos de 3 valores del blanco se elimina de la imagen
  #       copia[j,k,1:3] <- rep(NA, 3)
  #     }
  #   }
  # }
  #plot(copia)
  for(j in 1:3)
  {
    #se calculan las medianas requeridas solamente sobre los pixeles 
    #considerados como informativos de la propia imagen
    mediana <- median(copia[,,j], na.rm = TRUE)
    medianas[i, j] <- mediana 
  }
  remove(copia)
  remove(fruta)
}
#####se termina preprocesamiento
colnames(medianas) <- c('r', 'g', 'b')
#revise si habia diferencias en el encoding, porque la función da puntos en [0,1]
#mientras yo esperaba enteros en [0,255] y no hubo diferencias
reduce.a.medianas.discreta <- function(index)
{
  fruta <- load.image(file=imagenes[index])
  fruta[,,1:3] <- fruta[,,1:3]*255
  medianas <- apply(fruta[,,1:3], 3, median)
  medianas <- matrix(medianas, ncol = 1)
  row.names(medianas) <- c('r', 'g', 'b')
  return(medianas)
}
#se manipulan las medianas para convertirlas en un data.frame facil de manjer en las
#visualizaciones siguientes
imagenes.medianas <- medianas
df.imagenes.medianas <- as.data.frame(imagenes.medianas)
#no.normalizados <- lapply(1:length(imagenes), FUN= reduce.a.medianas.discreta)
#no.normalizados <- do.call("cbind", no.normalizados)
#df.no.normalizados <- as.data.frame(no.normalizados)
#df.no.normalizados <- as.data.frame(t(as.matrix(df.no.normalizados)))
#se obtienen los nombres de las imagenes provenientes del directorio
library(stringr)
imagenes.nombres <-str_split(imagenes,regex("[:digit:]"),simplify = T)[,1]# "Eliminar de un digito hacia adelante"
imagenes.nombres <- str_split(imagenes.nombres, regex("$"),simplify = T)[,1]  # Elimnar los "" finales
imagenes.nombres <-str_split(imagenes.nombres, regex("r$"),simplify = T)[,1]  # "Eliminar de un "_r" hacia adelante"
df.imagenes.medianas$tipo <- imagenes.nombres #se agregan al conjunto de datos los nombres 'tipos' de frutas
df.imagenes.medianas$tipo <- factor(df.imagenes.medianas$tipo)
df.imagenes.medianas$tipo2 <- df.imagenes.medianas$tipo
#en la siguiente seccion se crea una columna que agrupa por tipo de fruta sin importar su orientacion
levels(df.imagenes.medianas$tipo2)  <- c(rep('Apple_Braeburn',2), 
          'Apple_Golden', rep('Apple_Granny_Smith', 2),
          rep('Apricot', 2), rep('Avocado', 2), rep('Carambula', 2),
          rep('Cherry', 2), rep('Huckleberry', 2), rep('Kiwi', 2), rep('Orange', 2),
          rep('Peach', 2), rep('Pineapple', 2), rep('Strawberry', 2))

#se definen colores para las visualizaciones de los conjuntos
colores <- c('red', 'green', 'green4', 'pink4', 
             'black', 'yellow', 'purple', 'navy',
             'brown', 'orange', 'salmon', 
             'darkgoldenrod1', 'magenta')
#df.no.normalizados$tipo <- df.imagenes.medianas$tipo
#df.no.normalizados$tipo2 <- df.imagenes.medianas$tipo2
library(plotly) #visualizacion en 3D de las imagenes representadas por las medianas en cada canal
p1 <- plot_ly(df.imagenes.medianas, x = ~r, y = ~g,
             z = ~b, color = ~tipo2, 
             colors = colores) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'Red'),
                      yaxis = list(title = 'Green'),
                      zaxis = list(title = 'Blue')))
p1
#se realiza PCA sobre las medianas PCA con matriz de covarianzas
PCA <- princomp(df.imagenes.medianas[,1:3], cor = FALSE, scores = TRUE)
valores.pro <- PCA$sdev**2
cumsum(valores.pro)/sum(valores.pro) #las dos primeras componentes representan 96.73% de la varianza total
datos.rotados <- as.matrix(df.imagenes.medianas[,1:3])%*%PCA$loadings#se rotan los datos
datos.rotados <- as.data.frame(datos.rotados)
datos.rotados$tipo <- df.imagenes.medianas$tipo
datos.rotados$tipo2 <- df.imagenes.medianas$tipo2
datos.rotados$tipo3 <- colores[datos.rotados$tipo2]
p1 <- ggplot(datos.rotados, aes(x=Comp.1, y = Comp.2, 
                          color=tipo2)) +
  geom_point()+ 
  scale_color_manual(values=colores)+
  ggtitle('PCA imagenes') + theme_minimal() + 
  xlab('Primer componente principal') +
  ylab('Segunda componente principal') 
# se realiza kernel-pca
library(kernlab)
set.seed(0)
kernel.pca <- kpca(~., data=df.imagenes.medianas[,1:3], kernel="rbfdot", kpar=list(sigma=1/90))
datos.kernel <- as.data.frame(rotated(kernel.pca))
eigen.values <- eig(kernel.pca)
cumsum(eigen.values)/sum(eigen.values)#las dos primeras componentes explican el 96% de la varianz
datos.kernel$tipo2 <- df.imagenes.medianas$tipo2
p2 <- ggplot(datos.kernel, aes(x=V1, y = V2, 
                          color=tipo2)) +
  geom_point()+ 
  scale_color_manual(values=colores)+
  ggtitle('Kernel PCA (imagenes)') + theme_minimal() + 
  xlab('Kernel-Primer componente principal') +
  ylab('Kernel-Segunda componente principal') 
p2
library(ggpubr)
ggarrange(p1, p2, ncol = 2, nrow = 1)
#se realiza kmeans y kernel-kmeans
set.seed(0)
df.imagenes.medianas$kmeans <- factor(kmeans((df.imagenes.medianas[,1:3]), centers=length(colores), nstart = 100 )$cluster)
table(df.imagenes.medianas$tipo2, df.imagenes.medianas$kmeans)#error de 0.2461538
set.seed(0)
kk <- kkmeans(as.matrix((df.imagenes.medianas[,1:3])), centers = 13, 
        kernel='rbfdot',  kpar = list(sigma = 90))
df.imagenes.medianas$kkmeans <- factor(kk@.Data)
k1 <- ggplot(datos.rotados, aes(x=Comp.1, y = Comp.2, 
                         color=df.imagenes.medianas$kmeans)) +
  geom_point()+ 
  scale_color_manual(values=colores)+
  ggtitle('Grupos usando kmeans (imagenes)') + theme_minimal() + 
  xlab('Primer componente principal') +
  ylab('Segunda componente principal') 
kk1 <- ggplot(datos.rotados, aes(x=Comp.1, y = Comp.2, 
                          color=df.imagenes.medianas$kkmeans)) +
  geom_point()+ 
  scale_color_manual(values=colores)+
  ggtitle('Grupos usando kernel-kmeans (imagenes)') + theme_minimal() + 
  xlab('Primer componente principal') +
  ylab('Segunda componente principal') 
table( df.imagenes.medianas$tipo2, df.imagenes.medianas$kkmeans)
402/1300 #sigma de 100
473/1300 #sigma de 60
452/1300 #sigma de 200 
504/1300 #sigma de 80
388/1300 #sigma de 90
475/1300 #sigma de 95 
482/1300 #sigma de 85
408/1300 #sigma de 50
ggarrange(p1, k1, kk1, ncol = 2, nrow = 2)

############ inciso d imagenes.nombres
##de nuevo el preprocesamiento de quitar pixeles cercanos al blanco 
medianas.HSV <- matrix(rep(-1,9*index), nrow = index) #se reserva espacio para guardar las medianas de los incisos 1 a 3
for(i in 1:index)
{
  fruta <- load.image(file=imagenes[i])
  print(imagenes[i])
  copia <- fruta
  for(j in 1:100)#para cada pixel se mide su distancia con el pixel blanco
  {
    for(k in 1:100)
    {
      if(  L2(copia[j,k,1,1:3], rep(1.,3)) <=0.02037707)
      {
        #aproximadamente si el pixel dista en menos de 3 valores del blanco se elimina de la imagen
        copia[j,k,1, 1:3] <- rep(NA, 3)
      }
    }
  }
  #plot(copia)
  copia <- RGBtoHSV(copia)
  for(j in c(1,4,7))
  {
    #se calculan las medianas requeridas solamente sobre los pixeles 
    #considerados como informativos de la propia imagen
    canal <- round(j/3,0)+1
    mediana <- quantile(copia[,,1,canal], na.rm = TRUE, probs =c(.25,.5,.75)  )
    medianas.HSV[i, j:(j+2) ] <- mediana 
  }
  remove(copia)
  remove(fruta)
}
###### se termina preprocesamiento
#la siguinete funcion la utilice en su momento para extraer los
# cuartiles sin considerar el preprocesamiento dd quitar los pixeles cercanos al 
# pixel blanco, como se obtuvieron resultados pobres se abandono este esquema
# reduce.a.medianas.HSV <- function(index)
# {
#   fruta <- load.image(file=imagenes[index])
#   fruta <- RGBtoHSV(fruta)
#   mediana <- apply(fruta[,,1,1:3], 3, function(x){
#     quantile(x, probs = c(.25, .5, .75 ))
#   })
#   medianas <- matrix(mediana, byrow = TRUE, ncol = 9)
#   colnames(medianas) <- c('H.25','H.5', 'H.75',
#                           'S.25','S.5', 'S.75',
#                           'V.25','V.5', 'V.75')
#   row.names(medianas) <- index
#   return(medianas)
# }
imagenes.medianas.hsv <- medianas.HSV
#se manipulan los cuartiles para convertirlas en un data.frame facil de manjer en las
#visualizaciones siguientes
  #se copian columnas  on el tipo de fruta sin importar su orientacion y consideracndo su orientacion
df.imagenes.medianas.hsv <- as.data.frame(imagenes.medianas.hsv)
df.imagenes.medianas.hsv$tipo1 <- df.imagenes.medianas$tipo
df.imagenes.medianas.hsv$tipo2 <- df.imagenes.medianas$tipo2
#se realiza PCA sobre las medianas PCA con matriz de covarianzas
PCA <- princomp(df.imagenes.medianas.hsv[,1:9], cor = FALSE, scores = TRUE)
valores.pro <- PCA$sdev**2
cumsum(valores.pro)/sum(valores.pro) #las dos primeras componentes representan 99.9% de la varianza total
datos.rotados <- as.matrix(df.imagenes.medianas.hsv[,1:9])%*%PCA$loadings#se rotan los datos
datos.rotados <- as.data.frame(datos.rotados)
datos.rotados$tipo <- df.imagenes.medianas.hsv$tipo
datos.rotados$tipo2 <- df.imagenes.medianas.hsv$tipo2
p1 <- ggplot(datos.rotados, aes(x=Comp.1, y = Comp.2, 
                                color=tipo2)) +
  geom_point()+ 
  scale_color_manual(values=colores)+
  ggtitle('PCA imagenes') + theme_minimal() + 
  xlab('Primer componente principal') +
  ylab('Segunda componente principal') 
p1
# se realiza kernel-pca
library(kernlab)
set.seed(0)
kernel.pca <- kpca(~., data=df.imagenes.medianas.hsv[,1:9], 
                   kernel="rbfdot", kpar=list(sigma=.001*4))#.001*5
datos.kernel <- as.data.frame(rotated(kernel.pca))
eigen.values <- eig(kernel.pca)
cumsum(eigen.values)/sum(eigen.values)#las dos primeras componentes explican el 16.37% de la varianz
0.3597184
datos.kernel$tipo2 <- df.imagenes.medianas$tipo2
p2 <- ggplot(datos.kernel, aes(x=V1, y = V2, 
                               color=tipo2)) +
  geom_point()+ 
  scale_color_manual(values=colores)+
  ggtitle('Kernel PCA (imagenes)') + theme_minimal() + 
  xlab('Kernel-Primer componente principal') +
  ylab('Kernel-Segunda componente principal') 
p2
ggarrange(p1, p2, ncol = 2, nrow = 1)
#se realiza kmeans y kernel-kmeans
set.seed(0)
df.imagenes.medianas.hsv$kmeans <- factor(kmeans((df.imagenes.medianas.hsv[,1:9]), 
                                  centers=length(colores), nstart = 100 )$cluster)
table(df.imagenes.medianas.hsv$tipo2, df.imagenes.medianas.hsv$kmeans)#error de 0.2461538
401/1300 #error kmeans usando cuartiles
set.seed(0)
kk <- kkmeans(as.matrix(df.imagenes.medianas.hsv[,1:9]), 
              centers = 13, 
              kernel='rbfdot',  kpar = list(sigma = .04))
df.imagenes.medianas.hsv$kkmeans <- factor(kk@.Data)
k1 <- ggplot(datos.rotados, aes(x=Comp.1, y = Comp.2, 
                                color=df.imagenes.medianas.hsv$kmeans)) +
  geom_point()+ 
  scale_color_manual(values=colores)+
  ggtitle('Grupos usando kmeans (imagenes)') + theme_minimal() + 
  xlab('Primer componente principal') +
  ylab('Segunda componente principal') 
kk1 <- ggplot(datos.rotados, aes(x=Comp.1, y = Comp.2, 
                                 color=df.imagenes.medianas.hsv$kkmeans)) +
  geom_point()+ 
  scale_color_manual(values=colores)+
  ggtitle('Grupos usando kernel-kmeans (imagenes)') + theme_minimal() + 
  xlab('Primer componente principal') +
  ylab('Segunda componente principal') 
table( df.imagenes.medianas.hsv$tipo2, df.imagenes.medianas.hsv$kkmeans)
#con los cuartiles
#clasifica muy mal  #sigma de 100, 300,10
732/1300#sigma de 1
516/1300 #sigma de .1
475/1300 #sigma de .01
422/1300 #sigma de .05
434/1300 #sigma de .02
460/1300 #sigma de .03
422/1300 #sigma de .04
467/1300 #sigma de .06
540/1300 #sigma de .09
ggarrange(p1, k1, kk1, ncol = 2, nrow = 2)
#############contraste
colnames(df.imagenes.medianas.hsv) <- c('H.25', 'H.5', 'H.75',
                                    'S.25', 'S.5', 'S.75',
                                    'V.25', 'V.5', 'V.75',
                                    'tipo1', 'tipo2', 'kmeans', 'kkmeans')
library(plotly) #visualizacion en 3D de las imagenes representadas por las medianas en cada canal
p1 <- plot_ly(df.imagenes.medianas.hsv, x = ~H.75, y = ~S.75,
              z = ~V.75, color = ~tipo2, 
              colors = colores) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'H.75'),
                      yaxis = list(title = 'S.75'),
                      zaxis = list(title = 'V.75')))
p1

