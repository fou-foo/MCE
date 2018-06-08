####################################
#####    J. Antonio Garcia #########
#####   jose.ramirez@cimat.mx ######
####################################
## simulacion del dataset para ejemplificar el data pilling
setwd('/home/fou/Desktop/MCE/Second/CienciaDeDatos/DWD1')
n <- 20 
d <- 1000
library(MASS)
set.seed(0)
pos1 <- mvrnorm(n = n, mu =c(2.2, rep(0, d-1) ), Sigma = diag(rep(1,d)))
pos1 <- as.data.frame(pos1)
neg1 <- mvrnorm(n = n, mu =c(-2.2, rep(0, d-1) ), Sigma = diag(rep(1,d)))
neg1 <- as.data.frame(neg1)
save.image('data2.Rdata')
