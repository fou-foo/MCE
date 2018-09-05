# Tarea 1, Ejercicio 5.c y 5.d
# María Guadalupe Garrido Espinosa

#rm(list=ls())

min<-1
max<-10
n<- max-min+1 #Número de observaciones

x<- seq(min, max) #Espacio Muestral

samplesize<-10000 #El tamaño del que queremos la muestra

set.seed(13) #Fijamos la semilla en 13

y<-sample(x,size=samplesize,replace = TRUE, prob=NULL) 
#Generamos la muestra con reemplazo 
#(porque la muestra es de tamaño mayor al vector x)

freq<-table(y)#contamos 
results=as.data.frame(freq)#lo ponemos como una tabla de frecuencia

#Calculamos la media y la varianza
media=sum(y)/samplesize

varianza=sum( (y-media)^2)/samplesize

#Graficamos
plot(x, results$Freq, type = "h", col = "cornflowerblue", lwd = 3,ylim=c(0, max(results$Freq)),
     xlab="X",ylab="Frecuencia")


