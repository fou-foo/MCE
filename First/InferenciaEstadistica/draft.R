library(Matrix)
set.seed(0)
A <- Matrix(sample(1:16,16), 4, 4)
det(A)
p <- lu((A))
##
n <- dim(A)[1]
L <- diag(n)
L
A
#for ( k in 1:(n-1))
#{
  #for( i in (k+1):n )
  #{
    #L[i,k] <-  (A[i, k] / A[k, k] )
   # {
  #    for(j in (k):n )
     # {
    #    A[i, j] <- (A[i, j] - L[i,k]*A[k, j])
   #   }
  #  }
 # }
#}

#for(k in 1:(n-1))
#{
 # L[ (k+1):n,k] <- (A[(k+1):n, k]/A[k,k])
#  for(j in (k+1):n)
#  {
 #   A[j, ] <- (A[j, ] - L[j,k]*A[k, ] )
#  }
#}
expand(p)
det(L)
det(A)
A
for(k in 1:(n-1))
{
  L[ (k+1):n,k] <- (A[(k+1):n, k]/A[k,k])
  for(j in (k+1):n)
  {
    A[j, ] <- (A[j, ] - L[j,k]*A[k, ] )
  }
}
