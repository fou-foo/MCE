##############################################################
library(rpart)
mypred <- function(r1,...)
{
    #funcion para predecir recibe un objeto de la clase 'rpart' 
    #regresa la prediction maximo verosimil 
    y_hat <- predict(r1, ...)
    y_hat <- apply(y_hat, 1, which.max)
    y_hat
}
boost <- function(x, y, M, maxdepth=2, normalpha=F, cp = 0.01, minsplit = 50)
{
    #funcion de boost para arboles
    #x (dataframe)
    #y (vector) : corresponde a las etiquetas de las observaciones de 'x'
    N <- length(y)
    w <- rep(1/N, N) # apriori uniforme
    yhat <- matrix(0, ncol = M, nrow = N)
    err <-  rep(1, M)
    alpha <- rep(1/M, M) #pesos identicos al inicio
    treelist <- vector('list', M) #lista para cuardar los clasificadores debiles
    data.train <- data.frame(x=x,y=y)
    k <- length(unique(data.train$y)) #numero de categorias
    for(m in 1:M)
    {
        if(m%%50==0)
            print(m) #para darnos una idea de en que iteracion va 
        treelist[[m]] <- rpart(y~., data=data.train, method='class', 
                               control = rpart.control(cp = cp, minsplit=minsplit, maxdepth=maxdepth), 
                               weights=w) #entremos las funciones base
        yhat[,m] <- mypred(treelist[[m]]) #recalibramos
        diferentes <- yhat[,m]!= as.numeric(y)
        if(sum(diferentes) == 0) #por si en un conjunto facíl se logra el cero (como en el conjunto de datos 'iris')
        {
            indices <- lapply(treelist, FUN = function(x) {return(!is.null(x))})
            indices <- unlist(indices)
            cat('Error cero alcanzado en iteracion: ', m)
            return(list(treelist=treelist[indices], 
                        alpha=alpha[indices], 
                        weights=w[indices], error=err[indices]))
        }
        err[m] <- sum(w*( diferentes )) /sum(w) #numero de errores en este clasificador debil 'm'
        alpha[m] <- log((1-err[m])/err[m]) + log( k -1 )  #la GRAN diferencia entre arboles y adaboost
        w <- w*exp(alpha[m]* ( diferentes )) #actualizacion de pesos
        w <- w/sum(w**2)**.5
    }
    #if(normalpha)
    #   alpha <- alpha/sum(alpha)
    result <- list(treelist=treelist,alpha=alpha,weights=w,error=err)
    class(result) <- 'boost'
    return(invisible(result))
}
predict.boost <- function(bo, x, y, M=length(bo$treelist), verbose=F)
{
    newdata <- data.frame(x=x,y=y)
    yhat <- rep(-1,nrow(newdata))
    k <- length(unique(y))
    eval <- rep(-1, k)
    cuenta <- 0
    votos <- rep(-1,M)
    for(j in 1:nrow(newdata)) #para cada observacion estimamos su prediccion 
    {
        if(j%%100==0)
            print(j) #para darnos una idea de en que iteracion va 
        eval <- rep(-1, k)
        for(m in 1:M)#evaluamos la observacion en cada arbol
        {
            votos[m] <- mypred( bo$treelist[[m]], newdata=newdata[j,], weights=bo$weights ) 
        }
        for (l in 1:k) #guardamos los resultados de cada clase
        {
            eval[l] <- sum(bo$alpha*(votos == l))
        }
        yhat[j] <- which.max(eval) #regresamos la maximo verosimil
    }
    return(list( class=yhat))
}
############################################################
setwd('/home/fou/Desktop/MCE_CIMAT/Second/CienciaDeDatos/Tarea4')
###### uso de la implementacion en el dataset 'vowel'###########
set.seed(0)
data <- read.table('vowel.txt', header = FALSE)
colnames(data) <- c('Train.or.Test', 'Speaker.Number', 'Sex', 'Feature.0', 'Feature.1',
                    'Feature.2', 'Feature.3', 'Feature.4', 'Feature.5',
                    'Feature.6',  'Feature.7', 'Feature.8', 'Feature.9', 'Class')
apply(data, 2, class)
data$Speaker.Number <- factor(data$Speaker.Number)
data$Sex <- factor(data$Sex)
data$Class <- factor(data$Class)
train <- subset(data, Train.or.Test == 0)
test <- subset(data, Train.or.Test == 1)
y.train <- train$Class
train$Train.or.Test  <- train$Class <- NULL
y.test <- test$Class
test$Train.or.Test <-  test$Class <- NULL
set.seed(0)
adaboost <- boost(x = train, y = y.train , M = 400, maxdepth = 12, cp = -1, minsplit = 10)
#i600 <- data.frame(x = 1:600, error = adaboost$error)
#i600$iteraciones <- '600'
#i400 <- data.frame(x = 1:400, error =adaboost$error)
#i400$iteraciones <- '400'
#i200 <- data.frame(x = 1:200, error =adaboost$error)
#i200$iteraciones <- '200'
#errores.vowel <- rbind(i600, i400, i200)
#library(ggplot2)
#ggplot(errores.vowel, aes(x =x, y = error, color = iteraciones)) +geom_line() +
#    facet_grid(~iteraciones) + theme_minimal() + xlab('No. de calsificadores debiles') + ylab('Error en train')+
 #   ggtitle("Conjunto de datos 'Vowel'") + scale_color_manual(values=c('purple', 'orange', 'green4'))
plot(adaboost$error, type = 'l')
y_hat <- predict.boost(bo = adaboost, x = test, y = y.test)
sum(diag(table(y_hat$class, y.test)))/sum(table(y_hat$class, y.test))
#dos dos adaboost <- boost(x = train, y = y.train , M = 600, maxdepth = 10, cp = -1, minsplit = 10):#accuracy .5519481, 
#chido adaboost <- boost(x = train, y = y.train , M = 600, maxdepth = 12, cp = -1, minsplit = 10):#accuracy.5714286 
#esperando salida de lo que corre para compararla fijando semilla a cero derl modelo chido, con 200 en el train:0.5670996
# accuracy con 400 iteraciones: 554116
###################ejemplo rtificial  ###########################
n <- 100
u <- runif(300, 0, 1)
v1 <- max(6 - abs(1:100 - 11), 0)
v2 <- v1*(1:100-4)
v3 <- v1*(1:100+4)
x1 <- u[1:100]*v1+(1-u[1:100])*v2+rnorm(n)
x1 <- data.frame(x =x1, class='1')
x2 <- u[101:200]*v1 + (1-u[101:200])*v3+rnorm(n)
x2 <- data.frame(x=x2, class = '2')
x3 <- u[201:300]*v2+(1-u[201:300])*v3+rnorm(n)
x3 <- data.frame(x = x3, class = '3')
train <- rbind(x1, x2, x3)
y.train <- train$class
train$class <- NULL
set.seed(0)
n <- 1666
u <- runif(n*3, 0, 1)
v1 <- max(6 - abs(1:n - 11), 0)
v2 <- v1*(1:n-4)
v3 <- v1*(1:n+4)
x1 <- u[1:n]*v1+(1-u[1:n])*v2+rnorm(n)
x1 <- data.frame(x =x1, class='1')
x2 <- u[(n+1):2*n]*v1 + (1-u[(n+1):(2*n)])*v3+rnorm(n)
x2 <- data.frame(x=x2, class = '2')
x3 <- u[(2*n+1):(3*n)]*v2+(1-u[(2*n+1):(3*n)])*v3+rnorm(n)
x3 <- data.frame(x = x3, class = '3')
test <- rbind(x1, x2, x3)
y.test <- test$class
test$class <- NULL
set.seed(0)
adaboost <- boost(x = train, y = y.train , M = 600, maxdepth = 12, cp = -1, minsplit = 3)
plot(adaboost$error, type = 'l')
y_hat <- predict.boost(bo = adaboost, x = test, y = y.test)
sum(diag(table(y_hat$class, y.test)))/sum(table(y_hat$class, y.test))
i600 <- data.frame(x = 1:600, error = adaboost$error)
i600$iteraciones <- '600'
i400 <- data.frame(x = 1:400, error =adaboost$error)
i400$iteraciones <- '400'
i200 <- data.frame(x = 1:200, error =adaboost$error)
i200$iteraciones <- '200'
errores.vowel <- rbind(i600, i400, i200)
library(ggplot2)
ggplot(errores.vowel, aes(x =x, y = error, color = iteraciones)) +geom_line() +
    facet_grid(~iteraciones) + theme_minimal() + xlab('No. de calsificadores debiles') + ylab('Error en train')+
    ggtitle("Conjunto de datos 'Waveform'") + scale_color_manual(values=c('purple', 'orange', 'green4'))
#chafa adaboost <- boost(x = train, y = y.train , M = 600, maxdepth = 30, cp = -1, minsplit = 10)#accuracy:  
##chafa adaboost <- boost(x = train, y = y.train , M = 600, maxdepth = 20, cp = -1, minsplit = 9)#accuracy: .6786715
#se ve con mcadre adaboost <- boost(x = train, y = y.train , M = 600, maxdepth = 12, cp = -1, minsplit = 3)#accuracy (esperando este numero)
#accuracy con 200:0.3447379 validando 200 ok
#accuracy con 400: 0.3447 validando 
#accuracy con 600: 0.3447
#z <- data.frame(x = train, y = y.train)
#ggplot(z, aes(x=x, y = x, color=y))+geom_point()+facet_grid(~y)
plot(i600$)