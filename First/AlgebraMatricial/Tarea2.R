
menu <- function()
{
  print("Ingrese dimensión de la matriz")                   # Campuramos los datos para construir la matriz a factorizar
  n <- scan()
  print("Ingrese elementos de la matriz")
  ins <- scan()
  ins <- matrix(ins, n, n, byrow = TRUE)
  return(ins)
}

LU <- function(M)
{
  
  # Funcion para calcular la descomposión LU de la matriz 'M'
  # Entradas : 
              # M matriz de nxn
  # Salidas:
              # una lista cuyos elementos son
                  # L: una matriz triangular inferior de dimension nxn
                  # U: una matriz triangular superior de dimension nxn
  n <- dim(M)[1]                                          # Se obtiene la dimensión de la matriz de entrada
  L <- diag(n)                                            # se iniciliza una matriz identidad para guardar los pivotes
    for(k in 1:(n-1))                                     # el algoritmo itera sobre las columnas 
    {
      if(M[k,k] == 0)
      {
        stop("La descomposición LU encontro un cero, por lo que la descomposición no es posible")
      }
      if(abs(M[k, k]) < 1e-20 )
      {
        stop("La descomposición LU encontro un termino cercano a cero por lo que la descomposición puede ser poco exacta ")
      }
      L[ (k+1):n,k] <- (M[(k+1):n, k]/M[k,k])             # se calculan los pivotes
        for(j in (k+1):n)
        {
          M[j, ] <- (M[j, ] - L[j,k]*M[k, ] )             # se actualiza la matriz usando los pivotes
        }
    }
  return(list(L=L, U=M))
}

PLU <- function(M)
{
  # Funcion para calcular la descomposión LU de la matriz 'M'
  # Entradas : 
              # M matriz de nxn
  # Salidas:
              # una lista cuyos elementos son
                    # L: una matriz triangular inferior de dimension nxn
                    # U: una matriz triangular superior de dimension nxn
                    # P: una matriz de permutación de dimension nxn
  
  n <- dim(M)[1]                                        # se inicializa la dimensión de la matriz,
  indixes <- 1:n                                        # un vector para guardar las permutaciones efectuadas por el algoritmo
  P <- diag(n)                                          # una matriz para guardar la permutacion
  for (k in 1:(n-1))                                    # el algoritmo itera sobre las columnas 
  {
    if(M[k,k] == 0)
    {
      stop("La descomposición PLU encontro un cero, por lo que la descomposición no es posible")
    }
    if(abs(M[k, k]) < 1e-20 )
    {
      stop("La descomposición PLU encontro un termino cercano a cero por lo que la descomposición puede ser poco exacta ")
    }
    pivot.max <- which.max( abs( M[k:n, k] ) ) + (k-1)  # se localiza el elemento de mayor tamaño sobre las columnas para determinar el k-esimo pivote
    temp <- M[ k ,  ]                                   # se intercambian las columnas para que el pivote sea el maximo 
    M[k, ] <- M[pivot.max, ]
    M[pivot.max, ] <- temp
    temp.indice <- indixes[k]
    indixes[k] <- indixes[pivot.max]                    # se registra la permutación para construir la matriz de permutación
    indixes[pivot.max] <- temp.indice
    ij <- (k+1):n
    M[ij, k] <-  M[ij, k]/M[k, k]                       # se calculan los pivotes de la siguiente iteración
    M[ij, ij] <- M[ij, ij] - outer(M[ij, k], M[k, ij])  # se actualiza la matriz
  }
  L <- M                                                # Se prepara la salida de la función
  L[!lower.tri(M)] <- 0     
  L <- L + diag(n)
  M[lower.tri(M)] <- 0
  return(list(L = L, U = M, P=P[indixes,])) 
}

# Primer ejemplo de ejecución
M <- menu()
# Ejemplo de caso de prueba
# set.seed(50)
# M <- sample(1:100,100)
# M <- matrix(M, 10, 10, byrow = TRUE   )
(L.U. <- LU(M) )

# Verificamos que las matrices son quivalentes
det(M)
det(L.U.$U)
(P.L.U <- PLU(M))
# Verificamos que las matrices son quivalentes
(P.L.U$P %*% M ) - (round(P.L.U$L %*% P.L.U$U, 0))

# Segundo ejemplo de ejecución
M <- matrix(rep(1:5, 5), byrow = TRUE, nrow = 5)
(L.U. <- LU(M) )
(P.L.U <- PLU(M))
