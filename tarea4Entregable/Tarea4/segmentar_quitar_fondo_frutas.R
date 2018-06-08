setwd('/home/fou/Desktop/MCE_CIMAT/Second/CienciaDeDatos/Tarea4/Ejercicio3/frutas_propias/recortes')
imagenes <- dir()
index <- length(imagenes)
medianas <- matrix(rep(-1,3*index), nrow = index) #se reserva espacio para guardar las medianas de los incisos 1 a 3
