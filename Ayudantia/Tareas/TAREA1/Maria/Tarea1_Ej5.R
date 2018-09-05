# Tarea 1, Ejercicio 5.a
# María Guadalupe Garrido Espinosa

#rm(list=ls())


graf_unif<- function(a,b){
  min<-a
  max<-b
  n<- max-min+1 #Número de observaciones
  
  x<- seq(min, max) #Espacio Muestral
  
  fp<- rep(1/n,n) #obtenemos la probabilidad de cada uno de los valores del Espacio Muestral
  
  cdf<-(x-min+1)/n
  
  #Graficamos la función de masa 
  plot(x, fp, type = "h", col = "cornflowerblue", lwd = 3, main=paste("PMF Uniforme Discreta con n=",n))
  
  #plot(x, cdf,ylim=c(0, 1), type = "b", col = "cornflowerblue", lwd = 3, main=paste( "CDF Uniforme Discreta con n=",n) )
 
  #Graficamos la funcion de distribución acumulada
  plot(ecdf(x), col = "cornflowerblue", lwd = 3, main=paste( "CDF Uniforme Discreta con n=",n) )
  
  return(x) 
}

#Graficamos la distribución uniforme con n=5
graf_unif(1,5)

graf_unif(1,10)

graf_unif(1,50)


plot(x, cdf,ylim=c(0, 1), type = "b", col = "cornflowerblue", 
     add=FALSE,verticals=TRUE,
     lwd = 3, main=paste( "CDF Uniforme Discreta con n=",n) )

