# Tarea 1, Ejercicio 7
# María Guadalupe Garrido Espinosa

##################################
#Creamos una función que nos calcule el factorial
##################################

cal_fact<-function(numero){
  
  num_ent<-as.integer(numero)#Volvemos entero el número que ingresó el usuario
  
  if(num_ent<0)
  {
    print("Lo siento, ingresa un entero positivo")
  }
  else if (num_ent==0 | num_ent==1)
  {
    fact=1;
  } 
  else{
    fact=1;
    for (i in 1:num_ent)
    {
      fact=fact*i; 
    }
  }
  return(fact)
}

###################################
#Calculamos una función que nos diga la función de masa
###################################

pdf_bin<-function(n,p,x1,x2)
{
  acumulada<-0
  
  for (i in x1:x2)
  {
    value= cal_fact(n)/( cal_fact(i)*cal_fact(n-i) )*p^(i)*(1-p)^(n-i)
    #print(value)
    acumulada=acumulada+value
    #print(acumulada)
    
  }
  return(acumulada)
}


#Ahora calculamos lo que nos requiere el inciso a y b
#P(X=0)
pdf_bin(123,0.31,0,0)
dbinom(0,123,0.31,log=FALSE)

#P(X=123)
pdf_bin(123,0.31,123,123)
dbinom(123,123,0.31,log=FALSE)

#P(X=62)
pdf_bin(123,0.31,62,62)
dbinom(62,123,0.31,log=FALSE)

#P(0<=X<=10)
pdf_bin(123,0.31,0,10)
pbinom(10,123,0.31)-pbinom(0,123,0.31)



###################################
#Ahora vamos a hacer un código para calcular los cuantiles
###################################

cuan_bin<-function(size,prob,q)
{
  i<- 0
  p<- 0
  while(p<=q)
  {
   #print(size)
  #  print(prob)
   p=pdf_bin(size,prob,0,i)
   #print(p)
   i=i+1
  }
  return(i-1)
}

#Ahora obtenemos lo que el ejercicio nos pide
#¿Existe alguna función en R que calcule cuantiles? Sí, qbinom
#Calcular el cuantil de 0.25
cuan_bin(123,0.31,0.25)
qbinom(0.25,123,0.31)

#Calcular el cuantil de 0.5
cuan_bin(123,0.31,0.5)
qbinom(0.5,123,0.31)

#Calcular el cuantil de 0.75
cuan_bin(123,0.31,0.75)
qbinom(0.75,123,0.31)

