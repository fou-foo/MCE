#la primer parte es simularla desde la definicion y la comparo con la parametrización de R
set.seed(0)
alpha <- 1
beta <- 3
#n <- 100000
n <- 3000000
u <- runif(n, 0, 1)
Weibulls <- -beta * ( (log(1-u))^(1/alpha)  )
hist(Weibulls, freq = FALSE)
lines(density(Weibulls), col="green")
w_hat <- rweibull(n, alpha, beta)
hist(w_hat, freq = FALSE)
lines(density(w_hat), col="red")
          ##jala chido 
# QQ-plot
Weibulls.F_hat <- sort(Weibulls)  
Weibulls.F_hat <- (1:length(Weibulls.F_hat))/(length(Weibulls.F_hat))#cumsum(Weibulls.F_hat)/sum(Weibulls.F_hat)
plot(sort(Weibulls), Weibulls.F_hat, type = 'l', col = 'purple') #solo la acumulada Weibull
teoria <- pweibull(sort(Weibulls), alpha, beta)
plot(Weibulls.F_hat, teoria)
abline(a = 0, b = 1, col = 'blue')
#qqplot(Weibulls.F_hat, pweibull(Weibulls.F_hat, alpha, beta))
plot(Weibulls.F_hat, pnorm(sort(Weibulls), mean(Weibulls), sd(Weibulls)), ylim = c(0,1), xlim = c(0,1))
abline(a = 0, b = 1, col = 'blue')
qqnorm(Weibulls.F_hat )
abline(a = 0, b = 1, col = 'blue')


#z <- rweibull(n, alpha, beta)
#plot(z, pweibull(z, alpha, beta) )
#lines(sort(z), (1:length(z))/length(z), col = 'red')
#plot((1:length(z))/length(z), pweibull((1:length(z))/length(z), alpha, beta), col = 'red' , xlim = c(0,1), ylim = c(0,1))
#abline(0,1)

mediana <- sort(Weibulls)[floor(length(Weibulls)/2)]
Weibulls2 <- Weibulls
names(Weibulls2) <- Weibulls
moda <- dweibull(Weibulls, alpha, beta)[which.max(dweibull(Weibulls, alpha, beta))]
fun <- function(x)
{
  ((x[2] * (log(2))^(1/x[1])  - mediana)^2+ ( (x[2]*( (x[1]-1) / x[1] )^(1/x[1]) - moda   )^2) ) 
}
optim(par = c(3.,100), fn = fun)
