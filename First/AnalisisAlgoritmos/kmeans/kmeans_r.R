 
 

inicializacion <- function(k, n, columnas, columnas.vector, datos )
{       #  input: k-numero de centroides int
  #         n- numero de registros en el dataset
  #         columnas-numero de columnas con las que se trabaja,
  # y sus indices columnas.vector
  #         data set original
  # output: >inicializa los centroides
  muestra <- sample(1:n, k)
  centroides <- matrix(rep(0,k*columnas), ncol = columnas)
  for (i in 1:k)
  {
    centroides[i, ] <- unlist(datos[muestra[i], columnas.vector])
  }
  return(centroides)
}

l2 <- function(centroide, observacion)
{
  #Calculo individual de la distancia entre dos vectores
  # input: centroide-uno de los centroides
  #        observacion-uno de los registros del dataset
  return( sqrt( sum( (centroide-observacion) ^ 2 ) ) )
}

M.l2<- function(distancias, centroides, n ,k, data, columnas.vector )
{
  # Funcion para actualizar las distancias a los puntos
  #  input: k-numero de centroides int
  #         n- numero de registros en el dataset
  #  input: centroides
  #         data set original
  #         indices columnas.vector
  # output: distancias actualizadas 
  for( i in 1:n )
  {
    for(j in 1:k)
    {
      distancias[i,j] <- l2( data[i,columnas.vector], centroides[j, ] )
    }
  }
  return(distancias)
}

etiquetamiento <- function(distancias, n, k )
{
  # input : n - numero de observaciones
  #      : k - numero de observaciones
  # distancias : distancias actualizadas
  
  for( i in 1:n )
  {
    distancias[i, k+1] <- which.min(distancias[i, 1:k ])
  }
  return(distancias)
}


actualizacion.centroides <- function(data, distancias, k, centroides, columnas.vector )                #Funcion para mover los centroides 
  # Los actualiza y los devuelve
{
  for(i in 1:k)
  {
    centroides[i,] <- media(data, distancias, i, k, columnas.vector)
  }
  return(centroides)
}


media <- function(data, distancias, centro, k, columnas.vector)
{
  # Actualizacion de centroide en cada cluster
  
  subset <- data[ distancias[, k+1] == centro,   columnas.vector ]
  suma <- rep(0, length(columnas.vector))
  for(ii in 1:dim(subset)[1])
  {
    for(jj in 1:length(columnas.vector))
      suma[jj ]<- suma[jj]+ subset[ii,jj] 
  }
  return(suma/dim(subset)[1])
}


#seet up 
data <- iris                          #Datos con lo que se trabaja
set.seed(1)
##kernel
kmeans.foo <- function(data, klusters, iteraciones, columnas.vector= c(1,2,3,4))
{
  k <- klusters                                #Numero de kluster deseados 
  n <- dim(data)[1]                     #Numero de registros en el dataset
  #columnas.vector: columnas sin la de clasificacion 
  columnas <- dim(data)[2] - 1
  #inicializacion de centroides
  centroides <- inicializacion(k, n, columnas, columnas.vector, data )
  distancias <- matrix(rep(0, (k+1)*n), nrow  = n) #inicializacion de matriz de
  for( i in 1:iteraciones)
  {
    distancias <- M.l2(distancias, centroides, n ,k, data, columnas.vector) # particion de voronoi
    distancias <- etiquetamiento(distancias, n, k )         #asignacion de kluster
    centroides <- actualizacion.centroides(data, distancias, k, centroides, columnas.vector )
  }
  return(distancias[, k+1])
}



