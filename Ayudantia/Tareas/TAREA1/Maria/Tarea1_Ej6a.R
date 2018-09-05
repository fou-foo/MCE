# Tarea 1, Ejercicio 6
# María Guadalupe Garrido Espinosa

#rm(list=ls())

#res_moneda=c(1:2) #Generamos vector con las posibles salidas del lanzamiento de la moneda

#lanzamientos=sample(res_moneda,10,replace = TRUE) #Simulamos 10 veces este lanzamiento

#res_exp=data.frame(table(lanzamientos))#Hacemos el conteo

#res_exp[2,2]#Imprimimos el número de veces que se obtuvo águila (número 2)

#Los renglones 8,10 y 12 los podemos resumir en uno:
#data.frame(table(sample(res_moneda,10,replace = TRUE) ) )[2,2]


#############################################################
#¿Qué tan cargada está la moneda? 
# El valor será 0.5 si es el inciso a y 0.3 si es el inciso c
#############################################################

p=0.5; 
n=1000000;#¿Cuántas repeticiones quieres?

simular<-function(n){
  res_moneda=c(1:2)
  simulacion=rep(0,n)
  for (i in c(1:n)) 
    {
    simulacion[i]=data.frame(table(sample(res_moneda,10,replace = TRUE,prob =c(1-p,p)) ) )[2,2];
  
  }
  returnValue(simulacion);
}
simulacion=simular(n)#Repetimos el experimento 10^6

simulacion[1:3]#Mostramos los tres primeros resultados

resultados=data.frame(table(simulacion)) #contamos para obtener las frecuencias

#Graficamos las frecuencias de las águilas
plot(resultados$simulacion,resultados$Freq)

#Graficamos la probabilidad de que salga cierto número de águilas
plot(resultados$simulacion,resultados$Freq/n)

#Ahora vamos a graficar la función de masa de una B(10,0.5)

bin_res=matrix(0,11,2);
for(i in 0:10)
{
  bin_res[i,1]=i;
  bin_res[i,2]=dbinom(x=i,10,p);
}

rm(i)#Borramos el contador

res_binom=as.data.frame(bin_res)

colnames(res_binom)=c("Num_aguilas","Prop_bin")

colnames(resultados)=c("Num_aguilas","Prop_exp")

#creamos un arrglo que tenga los tres resultados con proporciones
comp=merge(res_binom,resultados, by="Num_aguilas", all=TRUE)

#Graficamos la función de masa de una B(10,p) encima del gráfico de proporciones
plot(comp$Num_aguilas,comp$Prop_exp/n,type="h",col="cornflowerblue",lwd = 10,xlab="# Águilas",ylab="Proporción")
par(bty = 'n')
points(comp$Num_aguilas,comp$Prop_bin,type="p",col="red",lwd = 10)

####
#Notas:
# 
#Inciso b) ¿Qué observa?
#
#       -Lo que se observa es que a medida que se incrementa el número de repeticiones
#       del experimento, el gráfico de proporciones tiende a ser más similar a la función
#       de masa de la distribución de una B(10,0.5). Adicionalmente, en la simulación
#       no se da el caso donde el número de águilas sea 0 o 10, esto se puede explicar 
#       por lo poco probable que son estos eventos, quizá si se incrementara el 
#       número de repeticiones se podría tener al menos una ocurrencia.

#       Además que se observa que el número de águilas que más se repite es el 5(=np)
#       
#
#Inciso c) ¿Qué observa?
#
#       -Al igual que en el inciso b, se observa es que a medida que se incrementa el 
#       número de repeticiones el experimento, el gráfico de proporciones tiende a ser
#       más similar a la función de masa de la distribución de una B(10,0.3).
#       -En la simulación no se da el caso que se tengan 0 o 10 águilas, al 
#       igual que en el inciso anterior, es posible que a medida que se incrementa el 
#       número de simulaciones se obtenga al menos una ocurrencia.
#   
#       Además que se observa que el número de águilas que más se repite es el 3 (=np)
####

