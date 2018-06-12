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
#save.image('data2.Rdata')
############ primera simulacion datos normales
set.seed(0)#semilla fuera de las funciones de simulacion 
simula.init <- function(n.mas, n.menos, n.test, d)
{
    library(MASS)
    #CLOSURE para simular obetner el errror de clasificación de MDP 
    #esta funcion regresa una funcion
    #fijamos los parametros y la muestra
    n.mas <- n.mas
    n.menos <- n.menos
    n.test <- n.test
    d <- d
    I <- diag(d)
    n <- n.mas + n.menos
    pos1 <- mvrnorm(n = n.mas, mu =c(2.2, rep(0, d-1) ), Sigma = diag(rep(1,d)))
    pos1 <- as.data.frame(pos1)
    neg1 <- mvrnorm(n = n.menos, mu =c(-2.2, rep(0, d-1) ), Sigma = diag(rep(1,d)))
    neg1 <- as.data.frame(neg1)
    train <- rbind(pos1,neg1)
    train$label <- 1
    train$label[(n.mas+1):(2*n.mas)] <- -1
    function(i)
    {
        library(MASS)
        #########genracion de conjunto test
        test.1 <- mvrnorm(n = n.test/2, mu =c(2.2, rep(0, d-1) ), Sigma = diag(rep(1,d)))
        test.1 <- as.data.frame(test.1)
        test.1$label <- 1
        test.2 <- mvrnorm(n = n.test/2, mu =c(-2.2, rep(0, d-1) ), Sigma = diag(rep(1,d)))
        test.2 <- as.data.frame(test.2)
        test.2$label <- -1
        test <- rbind(test.1, test.2)
        gc()
        ##############evaluacion de MDP 
        pos.mean <- apply(pos1, 2, mean)
        neg.mean <- apply(neg1, 2, mean)
        X <- ginv(cov(train[,  - (d+1)])) %*% (pos.mean - neg.mean)
        MDP <- X/sum(X**2)**.5 #normalizamos el vector 
        test.predic <- as.matrix(test[, 1:d])%*%MDP
        test$y_hat <-  ifelse(test.predic>=0, 1, -1)
        error.MDP <- 1-sum(diag(table(test$label, test$y_hat)))/n.test
        gc()
        ####################evaluacion de RLR
        library(glmnet)
        RLR <- glmnet(x = as.matrix(train[,1:d]), y = factor(train$label) , family = "binomial")
        #plot(RLR)
        y.RLR <- predict(RLR, newx = as.matrix(test[, 1:d]), type = "class", s =c(0.01)) #el paper dice que usaron este valor de lambda
        e <- table( y.RLR, test$label)
        error.RLR <- sum(diag(e))/n.test
        gc()
        ########################3 evaluacion de svm
        library(e1071)
        svm <- svm(factor(label) ~ ., data = train, kernel ='linear', cost = 100)
        y.hat <- predict(svm, test)
        error.svm <- 1 - sum(diag(table(y.hat, test$label)))/n.test
        gc()
        #################### evaluacion de DWD
        library(kerndwd)
        kern <- vanilladot()
        dwd <- kerndwd(x=train[, 1:d], train$label, kern, lambda = 100)
        y.hat <- predict.kerndwd(dwd,  kern=kern, newx = as.matrix(test[, 1:d]),
                                 x = as.matrix(train[, 1:d]), type='class')
        error.DWD <- 1 - sum(diag(table(y.hat, test$label)))/n.test
        gc()
        library(fastAdaboost)
        train$label <- factor(train$label)
        ada <- adaboost( label ~ . , data=train, 100)
        test$label <- factor(test$label)
        y.hat <- predict(ada, newdata = test, type ='class') 
        error.ada <- sum(diag(table(y.hat$class, test$label)))/n.test
        res <- c(error.MDP, i, d, error.RLR, error.svm, error.DWD, error.ada)
        names(res) <- c('MDP', 'corrida','dimension', 'RLR', 'SVM', 'DWD', 'Adaboost')
        gc()
        return(res)
    }
}
t1 <- Sys.time()
simula10 <- simula.init(n.mas = 25, n.menos=25, n.test =200, d=10 )
library(parallel)
simulacion.10 <- mcmapply( FUN=simula10, 1:100, mc.cores=4) 
simulacion.10 <- as.data.frame(t(simulacion.10))
t1 <- Sys.time() - t1
t1 #tiempo hasta 1600 solo para MDP con 100
setwd('/home/fou/Desktop/MCE/Second/CienciaDeDatos/DWD1')
saveRDS(simulacion.10, file='simula10.rds')
t1 <- Sys.time()
simula40 <- simula.init(n.mas = 25, n.menos=25, n.test =200, d=40 )
simulacion.40 <- mcmapply( FUN=simula40, 1:100, mc.cores = 4) 
simulacion.40 <- as.data.frame(t(simulacion.40))
saveRDS(simulacion.40, file='simula40.rds')
t1 <- Sys.time() - t1
t1 #tiempo hasta 1600 solo para MDP con 100
t1 <- Sys.time()
simula100 <- simula.init(n.mas = 25, n.menos=25, n.test =200, d=100 )
simulacion.100 <- mcmapply( FUN=simula100, 1:100, mc.cores=4 ) 
simulacion.100 <- as.data.frame(t(simulacion.100))
saveRDS(simulacion.100, file='simula100.rds')
t1 <- Sys.time() - t1
t1 #tiempo hasta 1600 solo para MDP con 100
t1 <- Sys.time()
simula400 <- simula.init(n.mas = 25, n.menos=25, n.test =200, d=400 )
simulacion.400 <- mcmapply( FUN=simula400, 1:100, mc.cores=4) 
simulacion.400 <- as.data.frame(t(simulacion.400))
saveRDS(simulacion.400, file='simula400.rds')
t1 <- Sys.time() - t1
t1 #tiempo  10, 40, 100, 400, 1600
t1 <- Sys.time()
simula1600 <- simula.init(n.mas = 25, n.menos=25, n.test =200, d=1600 )
simulacion.1600 <- mcmapply( FUN=simula1600, 1:100, mc.cores = 4 ) 
simulacion.1600 <- as.data.frame(t(simulacion.1600))
saveRDS(simulacion.1600, file='simula1600.rds')
t1 <- Sys.time() - t1
t1 #tiempo  10, 40, 100, 400, 1600
###############simulacion plot 1
simulacion <- rbind(simulacion.10, simulacion.40, simulacion.100, simulacion.400, simulacion.1600)
library(reshape2)
names(simulacion)
a <- melt(simulacion, id =c('corrida', 'dimension')) 
library(dplyr)
library(ggplot2)
names(a)
a %>% group_by( dimension, variable ) %>% 
    summarise(media =mean(value), sd =sd(value) ) -> a
p <- ggplot(a, aes(x=dimension, y=media, colour=variable)) + 
    geom_errorbar(aes(ymin=media-sd, ymax=media+sd), width=.1) +
    geom_line() + geom_point() + scale_x_continuous(trans='log10') + theme_minimal()+
    ylab('Error de test') + ggtitle('Multinormales (una coordenada diferente)')
p
#######################################
#######3########### 
########## simulacion2
set.seed(0)#semilla fuera de las funciones de simulacion 
simula.init2 <- function(n.mas, n.menos, n.test, d)
{
    library(MASS)
    #CLOSURE para simular obetner el errror de clasificación de MDP 
    #esta funcion regresa una funcion
    #fijamos los parametros y la muestra
    n.mas <- round(n.mas*.8)
    n.menos <- round(n.menos*.8)
    n.test <- round(n.test*.8)
    d <- d
    I <- diag(d)
    n <- n.mas + n.menos
    pos1 <- mvrnorm(n = n.mas, mu =c(2.2, rep(0, d-1) ), Sigma = diag(rep(1,d)))
    pos1 <- rbind(pos1, mvrnorm(n = round(n.mas*.1), mu =c(100, rep(0, d-1) ), Sigma = diag(rep(1,d))))
    pos1 <- rbind(pos1, mvrnorm(n = round(n.mas*.1), mu =c(500, rep(0, d-1) ), Sigma = diag(rep(1,d))))
    pos1 <- as.data.frame(pos1)
    neg1 <- mvrnorm(n = n.menos, mu =c(-2.2, rep(0, d-1) ), Sigma = diag(rep(1,d)))
    neg1 <- rbind(neg1, mvrnorm(n = round(n.mas*.1), mu =c(-100, rep(0, d-1) ), Sigma = diag(rep(1,d))))
    neg1 <- rbind(neg1, mvrnorm(n = round(n.mas*.1), mu =c(-500, rep(0, d-1) ), Sigma = diag(rep(1,d))))
    neg1 <- as.data.frame(neg1)
    train <- rbind(pos1,neg1)
    train$label <- 1
    train$label[(n.mas+1):(2*n.mas)] <- -1
    function(i)
    {
        library(MASS)
        #########genracion de conjunto test
        test.1 <- mvrnorm(n = round(.8*n.test/2), mu =c(2.2, rep(0, d-1) ), Sigma = diag(rep(1,d)))
        test.1 <- rbind(test.1, mvrnorm(n = round(.1*n.test/2), mu =c(100, rep(0, d-1) ), Sigma = diag(rep(1,d))))
        test.1 <- rbind(test.1, mvrnorm(n = round(.1*n.test/2), mu =c(500, rep(0, d-1) ), Sigma = diag(rep(1,d))))
        test.1 <- as.data.frame(test.1)
        test.1$label <- 1
        test.2 <- mvrnorm(n = n.test/2, mu =c(-2.2, rep(0, d-1) ), Sigma = diag(rep(1,d)))
        test.2 <- rbind(test.2, mvrnorm(n = round(.1*n.test/2), mu =c(-100, rep(0, d-1) ), Sigma = diag(rep(1,d))))
        test.2 <- rbind(test.2, mvrnorm(n = round(.1*n.test/2), mu =c(-500, rep(0, d-1) ), Sigma = diag(rep(1,d))))
        test.2 <- as.data.frame(test.2)
        test.2$label <- -1
        test <- rbind(test.1, test.2)
        ##############evaluacion de MDP 
        pos.mean <- apply(pos1, 2, mean)
        neg.mean <- apply(neg1, 2, mean)
        X <- ginv(cov(train[,  - (d+1)])) %*% (pos.mean - neg.mean)
        MDP <- X/sum(X**2)**.5 #normalizamos el vector 
        test.predic <- as.matrix(test[, 1:d])%*%MDP
        test$y_hat <-  ifelse(test.predic>=0, 1, -1)
        error.MDP <- 1-sum(diag(table(test$label, test$y_hat)))/n.test
        gc()
        ####################evaluacion de RLR
        library(glmnet)
        RLR <- glmnet(x = as.matrix(train[,1:d]), y = factor(train$label) , family = "binomial")
        #plot(RLR)
        y.RLR <- predict(RLR, newx = as.matrix(test[, 1:d]), type = "class", s =c(0.01)) #el paper dice que usaron este valor de lambda
        e <- table( y.RLR, test$label)
        error.RLR <- sum(diag(e))/n.test
        gc()
        ########################3 evaluacion de svm
        library(e1071)
        svm <- svm(factor(label) ~ ., data = train, kernel ='linear', cost = 100)
        y.hat <- predict(svm, test)
        error.svm <- 1 - sum(diag(table(y.hat, test$label)))/n.test
        gc()
        #################### evaluacion de DWD
        library(kerndwd)
        kern <- vanilladot()
        dwd <- kerndwd(x=train[, 1:d], train$label, kern, lambda = 100)
        y.hat <- predict.kerndwd(dwd,  kern=kern, newx = as.matrix(test[, 1:d]),
                                 x = as.matrix(train[, 1:d]), type='class')
        error.DWD <- 1 - sum(diag(table(y.hat, test$label)))/n.test
        gc()
        library(fastAdaboost)
        train$label <- factor(train$label)
        ada <- adaboost( label ~ . , data=train, 100)
        test$label <- factor(test$label)
        y.hat <- predict(ada, newdata = test, type ='class') 
        error.ada <- sum(diag(table(y.hat$class, test$label)))/n.test
        res <- c(error.MDP, i, d, error.RLR, error.svm, error.DWD, error.ada)
        names(res) <- c('MDP', 'corrida','dimension', 'RLR', 'SVM', 'DWD', 'Adaboost')
        gc()
        return(res)
    }
}
t1 <- Sys.time()
simula10 <- simula.init2(n.mas = 25, n.menos=25, n.test =200, d=10 )
library(parallel)
simulacion.10 <- mcmapply( FUN=simula10, 1:100, mc.cores=4) 
simulacion.10 <- as.data.frame(t(simulacion.10))
t1 <- Sys.time() - t1
t1 #tiempo hasta 1600 solo para MDP con 100
setwd('/home/fou/Desktop/MCE/Second/CienciaDeDatos/DWD1')
saveRDS(simulacion.10, file='simula2_10.rds')
t1 <- Sys.time()
simula40 <- simula.init2(n.mas = 25, n.menos=25, n.test =200, d=40 )
simulacion.40 <- mcmapply( FUN=simula40, 1:100, mc.cores = 4) 
simulacion.40 <- as.data.frame(t(simulacion.40))
saveRDS(simulacion.40, file='simula2_40.rds')
t1 <- Sys.time() - t1
t1 #tiempo hasta 1600 solo para MDP con 100
t1 <- Sys.time()
simula100 <- simula.init2(n.mas = 25, n.menos=25, n.test =200, d=100 )
simulacion.100 <- mcmapply( FUN=simula100, 1:100, mc.cores=4 ) 
simulacion.100 <- as.data.frame(t(simulacion.100))
saveRDS(simulacion.100, file='simula2_100.rds')
t1 <- Sys.time() - t1
t1 #tiempo hasta 1600 solo para MDP con 100
t1 <- Sys.time()
simula400 <- simula.init2(n.mas = 25, n.menos=25, n.test =200, d=400 )
simulacion.400 <- mcmapply( FUN=simula400, 1:100, mc.cores=4) 
simulacion.400 <- as.data.frame(t(simulacion.400))
saveRDS(simulacion.400, file='simula2_400.rds')
t1 <- Sys.time() - t1
t1 #tiempo  10, 40, 100, 400, 1600
t1 <- Sys.time()
simula1600 <- simula.init2(n.mas = 25, n.menos=25, n.test =200, d=1600 )
simulacion.1600 <- mcmapply( FUN=simula1600, 1:100, mc.cores = 4 ) 
simulacion.1600 <- as.data.frame(t(simulacion.1600))
saveRDS(simulacion.1600, file='simula2_1600.rds')
t1 <- Sys.time() - t1
t1 #tiempo  10, 40, 100, 400, 1600
###############simulacion plot 1
simulacion <- rbind(simulacion.10, simulacion.40, simulacion.100, simulacion.400, simulacion.1600)
a <- melt(simulacion, id =c('corrida', 'dimension')) 
a %>% group_by( dimension, variable ) %>% 
    summarise(media =mean(value), sd =sd(value) ) -> a
p <- ggplot(a, aes(x=dimension, y=media, colour=variable)) + 
    geom_errorbar(aes(ymin=media-sd, ymax=media+sd), width=.1) +
    geom_line() + geom_point() + scale_x_continuous(trans='log10') + theme_minimal()+
    ylab('Error de test') + ggtitle('Outlier mixture distribution')
p
##################################################################
######################### tercera simulacion
#################################################################
set.seed(0)#semilla fuera de las funciones de simulacion 
simula.init3 <- function(n.mas, n.menos, n.test, d)
{
    library(MASS)
    #CLOSURE para simular obetner el errror de clasificación de MDP 
    #esta funcion regresa una funcion
    #fijamos los parametros y la muestra
    n.mas <- n.mas
    n.menos <- n.menos
    n.test <- n.test
    d <- d
    I <- diag(d)
    n <- n.mas + n.menos
    pos1 <- mvrnorm(n = n.mas, mu =rep(50, d/2 ), Sigma = diag(rep(1,d/2)))
    pos1 <- cbind(pos1, pos1**(2))
    pos1 <- as.data.frame(pos1)
    neg1 <- mvrnorm(n = n.menos, mu =rep(50,  d/2) , Sigma = diag(rep(49,d/2)))
    neg1 <- cbind(neg1, neg1**(2)  )
    neg1 <- as.data.frame(neg1)
    train <- rbind(pos1,neg1)
    train$label <- 1
    train$label[(n.mas+1):(2*n.mas)] <- -1
    function(i)
    {
        library(MASS)
        #########genracion de conjunto test
        test.1 <- mvrnorm(n = n.test/2, mu =rep(50 , d/2 ), Sigma = diag(rep(1,d/2)))
        test.1 <- cbind(test.1, test.1**(2))
        test.1 <- as.data.frame(test.1)
        test.1$label <- 1
        test.2 <- mvrnorm(n = n.test/2, mu =rep(50, d/2) , Sigma = diag(rep(49,d/2)))
        test.2 <- cbind(test.2, test.2**(2))
        test.2 <- as.data.frame(test.2)
        test.2$label <- -1
        test <- rbind(test.1, test.2)
        ##############evaluacion de MDP 
        pos.mean <- apply(pos1, 2, mean)
        neg.mean <- apply(neg1, 2, mean)
        X <- ginv(cov(train[,  - (d+1)])) %*% (pos.mean - neg.mean)
        MDP <- X/sum(X**2)**.5 #normalizamos el vector 
        test.predic <- as.matrix(test[, 1:d])%*%MDP
        test$y_hat <-  ifelse(test.predic>=0, 1, -1)
        error.MDP <- 1-sum(diag(table(test$label, test$y_hat)))/n.test
        gc()
        ####################evaluacion de RLR
        library(glmnet)
        RLR <- glmnet(x = as.matrix(train[,1:d]), y = factor(train$label) , family = "binomial")
        #plot(RLR)
        y.RLR <- predict(RLR, newx = as.matrix(test[, 1:d]), type = "class", s =c(0.01)) #el paper dice que usaron este valor de lambda
        e <- table( y.RLR, test$label)
        error.RLR <- sum(diag(e))/n.test
        gc()
        ########################3 evaluacion de svm
        library(e1071)
        svm <- svm(factor(label) ~ ., data = train, kernel ='linear', cost = 100)
        y.hat <- predict(svm, test)
        error.svm <- 1 - sum(diag(table(y.hat, test$label)))/n.test
        gc()
        #################### evaluacion de DWD
        library(kerndwd)
        kern <- vanilladot()
        dwd <- kerndwd(x=train[, 1:d], train$label, kern, lambda = 100)
        y.hat <- predict.kerndwd(dwd,  kern=kern, newx = as.matrix(test[, 1:d]),
                                 x = as.matrix(train[, 1:d]), type='class')
        error.DWD <- 1 - sum(diag(table(y.hat, test$label)))/n.test
        gc()
        library(fastAdaboost)
        train$label <- factor(train$label)
        ada <- adaboost( label ~ . , data=train, 100)
        test$label <- factor(test$label)
        y.hat <- predict(ada, newdata = test, type ='class') 
        error.ada <- sum(diag(table(y.hat$class, test$label)))/n.test
        res <- c(error.MDP, i, d, error.RLR, error.svm, error.DWD, error.ada)
        names(res) <- c('MDP', 'corrida','dimension', 'RLR', 'SVM', 'DWD', 'Adaboost')
        gc()
        return(res)
    }
}
t1 <- Sys.time()
simula10 <- simula.init3(n.mas = 25, n.menos=25, n.test =200, d=10 )
library(parallel)
simulacion.10 <- mcmapply( FUN=simula10, 1:100, mc.cores=4) 
simulacion.10 <- as.data.frame(t(simulacion.10))
t1 <- Sys.time() - t1
t1 #tiempo hasta 1600 solo para MDP con 100
setwd('/home/fou/Desktop/MCE/Second/CienciaDeDatos/DWD1')
saveRDS(simulacion.10, file='simula3_10.rds')
t1 <- Sys.time()
simulacion.10 <- readRDS('simula3_10.rds')
simula40 <- simula.init3(n.mas = 25, n.menos=25, n.test =200, d=40 )
simulacion.40 <- mcmapply( FUN=simula40, 1:100, mc.cores = 4) 
simulacion.40 <- as.data.frame(t(simulacion.40))
saveRDS(simulacion.40, file='simula3_40.rds')
simulacion.40 <- readRDS('simula3_40.rds')
t1 <- Sys.time() - t1
t1 #tiempo hasta 1600 solo para MDP con 100
t1 <- Sys.time()
simula100 <- simula.init3(n.mas = 25, n.menos=25, n.test =200, d=100 )
simulacion.100 <- mcmapply( FUN=simula100, 1:100, mc.cores=4 ) 
simulacion.100 <- as.data.frame(t(simulacion.100))
saveRDS(simulacion.100, file='simula3_100.rds')
simulacion.100 <- readRDS('simula3_100.rds')
t1 <- Sys.time() - t1
t1 #tiempo hasta 1600 solo para MDP con 100
t1 <- Sys.time()
simula400 <- simula.init3(n.mas = 25, n.menos=25, n.test =200, d=400 )
simulacion.400 <- mcmapply( FUN=simula400, 1:100, mc.cores=4) 
simulacion.400 <- as.data.frame(t(simulacion.400))
saveRDS(simulacion.400, file='simula3_400.rds')
simulacion.400 <- readRDS('simula3_400.rds')
t1 <- Sys.time() - t1
t1 #tiempo  10, 40, 100, 400, 1600
t1 <- Sys.time()
simula1600 <- simula.init3(n.mas = 25, n.menos=25, n.test =200, d=1600 )
simulacion.1600 <- mcmapply( FUN=simula1600, 1:100, mc.cores = 4 ) 
simulacion.1600 <- as.data.frame(t(simulacion.1600))
saveRDS(simulacion.1600, file='simula3_1600.rds')
simulacion.1600 <- readRDS('simula3_1600.rds')
t1 <- Sys.time() - t1
t1 #tiempo  10, 40, 100, 400, 1600
###############simulacion plot 1
simulacion <- rbind(simulacion.10, simulacion.40, simulacion.100, simulacion.400, simulacion.1600)
a <- melt(simulacion, id =c('corrida', 'dimension')) 
a %>% group_by( dimension, variable ) %>% 
    summarise(media =mean(value), sd =sd(value) ) -> a
p <- ggplot(a, aes(x=dimension, y=media, colour=variable)) + 
    geom_errorbar(aes(ymin=media-sd, ymax=media+sd), width=.1) +
    geom_line() + geom_point() + scale_x_continuous(trans='log10') + theme_minimal()+
    ylab('Error de test') + ggtitle('Distribución de esferas anidadas')
p

