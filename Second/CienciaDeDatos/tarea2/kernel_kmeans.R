######### Implementacion de kernel k-means con shift basado en el paper:
#########Inderjit Dhillon, Yuqiang Guan and Brian Kulis.
#####A Unified view of Kernel k-means, Spectral Clustering and Graph Cuts.
Kernel.normal <- function(x,y, 
                         # alpha = (.05/2)**.5) #kernel gaussiano
                          alpha = 0.0013984457316839 ) #para el caso de vinos
{
  return(exp(-sum((x-y)**2)/(2*alpha**2)))
}
Kernel.poli <- function(x,y,c = 1, d =2) #kernel polinomial
{
  return((sum(x*y)+c)**d)
}
Kernel.kmeans.init <- function(data, cols,  f=Kernel.normal )
{
  # data (data.frame): data set con las observaciones
  # cols (vector numeric): indices de las columnas para el calculo 
  #de las distancias de 'data'
  #f (function): kernel a usar
  #se genera un dataframe con los pares de indices de las observaciones
  #con la finalidad de agilizar el calculo de la matriz de distancias 'MK'
  # ESTA FUNCION ES UN CLOSURE, REGRESA OTRA FUNCION, funciona como un constructor de clase del paradigma POO
  # LA FINALIDAD ES CALCULUAR LA MATRIZ DE KERNEL UNA SOLA VEZ Y PROBAR DIFERENTES VALORES DE PARAMETROS
  n <- dim(data)[1]
  indices <- 1:n
  #modificacion del 31 de marzo solo calcular la diagonal superior de la matriz K
  index <- data.frame(x1 = rep(indices,times=n:1))
  x2 <-lapply(FUN = function(x){
    return(n:x)
  }, 1:n)
  x2 <- unlist(x2)
  index$x2 <-x2
  index <- index[ order(index$x1, index$x2),]
  library(parallel)#para incrementar la velocidad usaremos calculo multicore
  #comienza calculo de la matriz superior del kernel entre todos los pares de observaciones
  kernel.sup <- lapply(FUN=function(i, cols){
    f(x=data[index[i,1],][cols],
      y=data[index[i,2],][cols])
  }, X = 1:dim(index)[1], cols)#,  mc.cores = detectCores()-2 )
  kernel.sup <- unlist(kernel.sup)
  MK <- matrix(rep(0, n*n), ncol = n)
  for(i in 1:dim(index)[1])
  {
    x <- index[i, 1]
    y <- index[i, 2]
    MK[x, y] <- kernel.sup[i]
    MK[y, x] <- kernel.sup[i]
  }
  #termina calculo de matriz de kernel
  function(data, sigma,t, k )
  {
    #data (data.frame) con las observaciones a clasificar
    #sigma (numeric): shift mencionado en el paper
    #t (numeric): numero de iteraciones 
    #k (numeric): numero de clusters
    data$cluster <- sample(1:k,dim(data)[1], replace = TRUE)#asignacion inicial aleatoria
    for(x in 1:t)
    {
      #siguiendo el paper citado realizamos las iteraciones
      #hacemos uso del shift 
      nuevo.cluster <- lapply(FUN=function(i)
      {
        #parte del calculo que no depende del shift
        d <- rep(-1, k)
        for(z in 1:k)
        {
          j <- which(data$cluster == z)
          d[z]<- MK[i,i] - 2*sum(MK[i, j])/length(j) + 
            sum(MK[j,j])/(length(j))**2
          d[is.na(d)] <- 0 #para cachar los errores de underflow
        }
        pis <- data.frame(table(data$cluster))
        # suma del shift a las entradas correspondientes
        for(z in 1:k)
        {
          if(z !=data[i,'cluster'])
          {
            d[z] <- d[z] + sigma +sigma/pis[z,'Freq'] 
          }else{
            d[z] <- d[z] + sigma -sigma/pis[data[i,'cluster'],'Freq'] 
          }
        }
        d[is.na(d)] <- 0
        mini <- min(d)
        candidatos <- which(d==mini)
        c.asig <- sample(candidatos, 1)#romper empates aleatoriamente
        (1:K)[c.asig]#nueva asignacion de cluster para la observacion 
      }, X =1:dim(data)[1]#, mc.cores = detectCores()-2 
      )
      data$cluster <- unlist(nuevo.cluster) #juntamos los resultados
    }
    return(data$cluster) #cluster finales
  }
}
############################################################
set.seed(0)
#####simulacion de datos parecidos a los del paper mencionado
#####son dos circunferencias con centro (.5,.5) y radios 1 y 4
#####se agrega en cada eje ruido ~ N(0,sigma=1/10 ) y N(0,1/10)
r <- 1 #radio
n <- 100 #la cuarta parte del numero de puntos que se van a generar
#se genera la primer circunferencia con ruido
x <- seq(-r, r, length=n) 
y1 <- sqrt(r**2-x**2) + rnorm(n,0,r/10)
y2 <- -sqrt(r**2-x**2) - rnorm(n,0,r/10)
m.a1 <- data.frame(x=rep(x+.5, 2), y = c(y1+.5,y2+.5), clase=1)
#se genera la segunda circunferencia con ruido
r <- 4
x <- seq(-r, r, length=n)
y1 <- sqrt(r**2-x**2) + rnorm(n,0,r/40)
y2 <- -sqrt(r**2-x**2) - rnorm(n,0,r/40)
m.a2 <- data.frame(x=rep(x+.5, 2), y = c(y1+.5,y2+.5), clase=2)
m.a <- rbind(m.a1, m.a2) #nuestro primer conjunto de prueba
library(ggplot2)
ggplot(m.a,#visualizamos nuestro primer conjunto de prueba 
       aes(x=x, y=y, color = factor(clase))) + geom_point() +
    theme_minimal() + theme(legend.position='none') +
    ggtitle('Muestra aleatoria generada (400 obs)') +xlab('') + ylab('')
set.seed(0)
label<- kmeans(m.a, centers = 2, nstart = 100) #comparamos el desempeño de kmeans en vista 
#                              #de que apriori sabemos que son 2 grupos 
p1 <- ggplot(m.a,#visualizamos nuestro primer conjunto de prueba 
       aes(x=x, y=y, color = factor(label$cluster))) + geom_point() +
    theme_minimal() + theme(legend.position='none') +
    ggtitle('Agrupamiento de kmeans (accuracy .49%)') +xlab('') + ylab('')
#visualizamos los resultados 
              #de clasificacion usando kmeans
sum(diag(as.matrix(table(m.a$clase, label$cluster))))/dim(m.a)[1] #se calcula el accuracy que es de .265, .735,.4975
#uso sobre 
set.seed(0)
Kernel.kmeans.simu <- Kernel.kmeans.init(m.a, cols=1:2) 
labels <- Kernel.kmeans.simu(data= m.a, sigma = -1, t = 20, k =2)
m.a$cluster <- labels
p3 <- ggplot(m.a,#visualizamos nuestro primer conjunto de prueba 
             aes(x=x, y=y, color = factor(cluster))) + geom_point() +
    theme_minimal() + theme(legend.position='none') +
    ggtitle('Agrupamiento de kernel.kmeans (accuracy .72%)') +xlab('') + ylab('')
sum(diag(as.matrix(table(m.a$clase,m.a$cluster))))/dim(m.a)[1]
# accuracy .3275 con shif de -1, .2175 con -.05, .47 con -3, .4875 con -5 totos los casos con t = 10 , .68 con 0 
#t=20: .69 con -1,.7175 con  -.05, -.5575 con -3, -.4525 con -5, .2575 con 0
#t=100, .72  con-1,.2525 con -.05,.5125 con -3,.4925 con -5, .245 con 0
####################
#########Comparacion de resultados usando la implementacion de kernlab
library(kernlab)
set.seed(10)
kkmeans.m.a <- kkmeans(as.matrix(m.a[,1:2]), centers=2 , 
                       kernel = "rbfdot", 
                       kpar = list(sigma = (.05/2)**.5),
                       alg="kkmeans")#mismo kernel y parametro del kernel 
sum(diag(as.matrix(table(m.a$clase, kkmeans.m.a@.Data))))/dim(m.a)[1]#accuracy de 1
p2 <- ggplot(m.a, aes(x=x, y=y, color = factor(kkmeans.m.a@.Data)))+
    geom_point() +theme_minimal()+theme(legend.position='none') +
    ggtitle('Agrupamiento perfecto de kkmeans')+xlab('')+ylab('')
library(ggpubr)
ggarrange(p1,p2,p3, nrow = 3)
##########
setwd('C:\\Users\\fou-f\\Desktop\\MCE\\Second\\CienciaDeDatos\\tarea1')
vinos <- read.csv('wine_quality.csv')
set.seed(0)
muestra.vinos <- sample(1:dim(vinos)[1], 1300)
vinos <- vinos[muestra.vinos,]
vinos$kmeans <- kmeans(vinos[,2:12], centers=7, nstart = 100)$cluster
(table(vinos$quality, vinos$kmeans))
276/1300#.26
pca.vinos <- princomp(vinos[,2:12], cor = TRUE)
pca.vinos.prin <- pca.vinos$scores[,1:2]
pca.vinos.prin <- as.data.frame(pca.vinos.prin)
p1 <- ggplot(pca.vinos.prin,#visualizamos nuestro primer conjunto de prueba 
             aes(x=Comp.1, y=Comp.2, color = factor(vinos$kmeans))) + geom_point() +
    theme_minimal() + theme(legend.position='none') +
    ggtitle('Agrupamiento de kmeans (accuracy 21%)') +
    xlab('Primer componente principal') + ylab('Segunda componente principal')
#
set.seed(0)
t1 <- Sys.time()
Kernel.kmeans.vinos <- Kernel.kmeans.init(vinos, cols=2:12)
t1 <- Sys.time()-t1
t2 <- Sys.time()
labels <- Kernel.kmeans.vinos(data=vinos, sigma = -1, t=20, k=7)
t2 <- Sys.time() - t2
vinos$kernelmio <- labels
table(vinos$quality,vinos$kernelmio) #acurracy con shiff -1 .16
550/1300 #t = 20
548/1300 #t=100
p3 <- ggplot(pca.vinos.prin,#visualizamos nuestro primer conjunto de prueba 
       aes(x=Comp.1, y=Comp.2, color = factor(vinos$kernelmio))) + geom_point() +
    theme_minimal() + theme(legend.position='none') +
    ggtitle('Agrupamiento de Kernel.kmeans (accuracy 42%)') +
    xlab('Primer componente principal') + ylab('Segunda componente principal')
####################################################
##############comparación para el dataset de vinos
set.seed(0)
kkmeans.vinos <- kkmeans(x=as.matrix(vinos[,2:12]), centers=7 ,
                         kernel = "rbfdot",
                         alg="kkmeans")
table(vinos$quality, kkmeans.vinos@.Data)
302/1300
#.12 con el sigma que yo doy .1078, .1237 con el default que vale 0.00153861608749537 
#visualizamos el agrupamiento en dos dimensiones
p2 <- ggplot(pca.vinos.prin, aes(x=Comp.1, y=Comp.2,
                           color = factor(kkmeans.vinos@.Data)))+
    geom_point()+theme_minimal()+ theme(legend.position='none') +
    ggtitle('Agrupamiento de kkmeans (accuracy 23%)') +
    xlab('Primer componente principal') + ylab('Segunda componente principal')
ggarrange(p1,p2,p3, ncol= 3)
