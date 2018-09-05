# Tarea 1, Ejercicio 8
# María Guadalupe Garrido Espinosa

#rm(list=ls())

# Nota, este script tiene un tiempo de ejecución entre 7 y 8 minutos

n=1000000;#¿Cuántas repeticiones quieres?

simular<-function(n){
  urna=c(rep(1,46),rep(2,49))#El uno son bolas grises, el 2 bolas blancas
  simulacion=rep(0,n)
  for (i in c(1:n)) 
  {
    simulacion[i]=data.frame(table(sample(urna,20,replace = FALSE) ) )[1,2];
    
  }
  returnValue(simulacion);
}
simulacion=simular(n)#Repetimos el experimento 10^6

simulacion[1:3]#Mostramos los tres primeros resultados

resultados=data.frame(table(simulacion)) #contamos las frecuencias

#Graficamos las frecuencias de las bolas grises obtenidas
plot(resultados$simulacion,resultados$Freq)

#Graficamos la probabilidad de que salga cierto número de bolas grises
plot(resultados$simulacion,resultados$Freq/n)

#Ahora vamos a graficar la función de masa de una distribución hipergeométrica

hip_res=matrix(0,21,2);
for(i in 0:20)
{
  hip_res[i,1]=i;
  hip_res[i,2]=dhyper(x=i,46,49,20);
}
rm(i)#Borramos el contador

res_hip=as.data.frame(hip_res)

colnames(res_hip)=c("Num_grises","Prop_hip")

colnames(resultados)=c("Num_grises","Prop_exp")

#creamos un arrglo que tenga los tres resultados con proporciones
comp=merge(res_hip,resultados, by="Num_grises", all=TRUE)

#¿Cuál es la probabilidad de que al extraer 20 bolas de la urna 5 sean grises?
#Si lo calculamos con la simulación es:
comp[which(comp[,1]==5),3]/n

#Si lo calculamos con el cálculo de la binomial:
comp[which(comp[,1]==5),2]

#Graficamos uno encima del otro
plot(comp$Num_grises,comp$Prop_exp/n,type="h",col="cornflowerblue",lwd = 10,xlab="# Bolas Grises",ylab="Proporción")
par(bty = 'n')
points(comp$Num_grises,comp$Prop_hip,type="p",col="red",lwd = 10)



