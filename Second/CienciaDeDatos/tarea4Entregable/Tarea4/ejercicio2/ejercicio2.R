###Lectura de datos y formación d eun solo dataset para hacer cv: DATOS de DIGITOS 
setwd('/home/fou/Desktop/MCE_CIMAT/Second/CienciaDeDatos/Tarea4/ejercicio2/')
test <- read.table('oef.test')
train <- read.table('oef.train')
data <- rbind(test,train)
#saveRDS(data, 'data.rds')
##############################################################33
library(caret)
library(plotly)
library(ggplot2)
library(ripa)
library(parallel)
set.seed(0)
################ control para las 10 k-folds ###########################33
###################33 eigenfaces ########################################
# de la tarea 1 dije que usar 50 componentes era suficiente
train.lm.cv <- function(cv=10, data)
{
    #funcion para evaluar por cv lm
    set.seed(0)
    #data <- data[ sample(1:dim(data)[1], dim(data)[1]),]
    p <- 50
    bloque <- round(dim(data)[1]/10)
    breaks <- c(seq(1, dim(data)[1], by = bloque  ), dim(data)[1])
    y <- data[,1]
    y <- factor(y)
    Y <- model.matrix(~y-1)
    data$V1 <- NULL
    acc.train <- rep(0, cv)
    acc.test <- rep(0, cv)
    for (i in 1:(length(breaks)-1))
    {
        indices <- breaks[i]:breaks[i+1]
        test <- data[ indices,   ]
        train <- data[ -indices,   ]
        pca <- princomp(train) #como los datos ya estan escalados en [-1, 1] uso la matriz de varianzas y covarianzas
        z <- pca$loadings[,1:p]
        Z <- as.matrix(train) %*%z
        y_hat <- lm(Y[ -indices,] ~ . -1 , data = as.data.frame(Z))
        b <- y_hat$coefficients
        y_hat <- Z%*%b
        y.ouput <-  apply(y_hat, 1, which.max)
        y.ouput <- y.ouput-1
        Matrix.C.train <- caret::confusionMatrix( factor(y.ouput), y[-indices]) #error de train
        Y.test.hat <- (as.matrix(data[indices,])%*%z)%*%b
        res2 <- apply(Y.test.hat, 1, which.max)
        res2 <- res2 - 1
        Matrix.C.test <- caret::confusionMatrix( factor(res2), y[indices])
        acc.train[i] <- Matrix.C.train$overall['Accuracy']
        acc.test[i] <-  Matrix.C.test$overall['Accuracy']
    }
    return( list(train=acc.train, test =acc.test) )
}
error.lm <- train.lm.cv(cv = 10, data)
error.lm  <- as.data.frame(error.lm)
error.lm$modelo <- 'lm'
##############vamos por el lda#############################################
train.LD_da.cv <- function(cv=10, data)
{
    library(MASS)
    #funcion para evaluar por cv lm
    set.seed(0)
    bloque <- round(dim(data)[1]/cv)
    breaks <- c(seq(1, dim(data)[1], by = bloque  ), dim(data)[1])
    acc.train <- rep(0, cv)
    acc.test <- rep(0, cv)
    for (i in 1:(length(breaks)-1))
    {
        indices <- breaks[i]:breaks[i+1]
        test <- data[ indices,   ]
        train <- data[ -indices,   ]
        y_hat <- lda( factor(V1) ~ ., data = data[-indices,]) 
        y.ouput <-  predict(y_hat, data[-indices, ])
        y.ouput <- y.ouput$class
        Matrix.C.train <- caret::confusionMatrix( factor(y.ouput), factor(data$V1[-indices])) #error de train
        Y.test.hat <- predict(y_hat, data[indices,])
        res2 <- Y.test.hat$class
        Matrix.C.test <- caret::confusionMatrix( factor(res2), factor(data$V1[indices]))
        acc.train[i] <- Matrix.C.train$overall['Accuracy']
        acc.test[i] <-  Matrix.C.test$overall['Accuracy']
    }
    return( list(train=acc.train, test =acc.test) )
}
error.lda <- train.LD_da.cv(cv = 10, data)
error.lda <- as.data.frame(error.lda)
error.lda$modelo <- 'lda'
###############vamos por el QDA   ########################################
train.QDA.cv <- function(cv = 10, data)
{
    #funcion para evaluar por cv qda
    set.seed(0)
    #data <- data[ sample(1:dim(data)[1], dim(data)[1]),]
    p <- 50
    bloque <- round(dim(data)[1]/10)
    breaks <- c(seq(1, dim(data)[1], by = bloque  ), dim(data)[1])
    y <- data[,1]
    y <- factor(y)
    data$V1 <- NULL
    acc.train <- rep(0, cv)
    acc.test <- rep(0, cv)
    for (i in 1:(length(breaks)-1))
    {
        indices <- breaks[i]:breaks[i+1]
        test <- data[ indices,   ]
        train <- data[ -indices,   ]
        pca <- princomp(train) #como los datos ya estan escalados en [-1, 1] uso la matriz de varianzas y covarianzas
        z <- pca$loadings[,1:p]
        Z <- as.matrix(train) %*%z
        y_hat <- qda(y[ -indices] ~ .  , data = as.data.frame(Z))
        y.ouput <-  predict(y_hat, as.data.frame(Z))
        y.ouput <- y.ouput$class
        Matrix.C.train <- caret::confusionMatrix( factor(y.ouput), factor(y[-indices])) #error de train
        Z.y <- as.matrix(test) %*% z
        Y.test.hat <- predict(y_hat, as.data.frame(Z.y))
        res2 <- Y.test.hat$class
        Matrix.C.test <- caret::confusionMatrix( factor(res2), factor(y[indices]))
        acc.train[i] <- Matrix.C.train$overall['Accuracy']
        acc.test[i] <-  Matrix.C.test$overall['Accuracy']
    }
    gc()
    return( list(train=acc.train, test =acc.test) )
}
error.qda <- train.QDA.cv(cv = 10, data)
error.qda <- as.data.frame(error.qda)
error.qda$modelo <- 'qda'
################ vamos por la red ###############################3
train.NN <- function(cv=10, data)
{
    #funcion para evaluar por cv NN
    library(nnet)
    set.seed(0)
    bloque <- round(dim(data)[1]/cv)
    breaks <- c(seq(1, dim(data)[1], by = bloque  ), dim(data)[1])
    acc.train <- rep(0, cv)
    acc.test <- rep(0, cv)
    for (i in 1:(length(breaks)-1))
    {
        indices <- breaks[i]:breaks[i+1]
        test <- data[ indices,   ]
        train <- data[ -indices,   ]
        y_hat <- nnet(factor(V1)  ~ ., data = data[-indices,], size = 3) #siguiendo la recomendación de hestie
        y.ouput <-  predict(y_hat, data[-indices, ], type = 'class')
        Matrix.C.train <- caret::confusionMatrix( factor(y.ouput), factor(data$V1[-indices])) #error de train
        Y.test.hat <- predict(y_hat, data[indices,], type = 'class')
        res2 <- Y.test.hat
        Matrix.C.test <- caret::confusionMatrix( factor(res2), factor(data$V1[indices]))
        acc.train[i] <- Matrix.C.train$overall['Accuracy']
        acc.test[i] <-  Matrix.C.test$overall['Accuracy']
    }
    gc()
    return( list(train=acc.train, test =acc.test) )
}
error.nn <- train.NN(cv = 10, data)
error.nn <- as.data.frame(error.nn)
error.nn$modelo <- 'nn'
################ vamor por el favorito de todos SVM###########################
train.SVM.init <- function(cv = 10, data, x)
{
    library(e1071)
    set.seed(0)
    bloque <- round(dim(data)[1]/cv)
    breaks <- c(seq(1, dim(data)[1], by = bloque  ), dim(data)[1])
    acc.train <- rep(0, cv)
    acc.test <- rep(0, cv)
    function(x){
        for (i in 1:(length(breaks)-1))
        {
            indices <- breaks[i]:breaks[i+1]
            test <- data[ indices,   ]
            train <- data[ -indices,   ]
            y_hat <- svm(factor(V1)  ~ ., data = data[-indices,], kernel='linear', cost = x ) #siguiendo la recomendación de hestie
            y.ouput <-  predict(y_hat, data[-indices, ], type = 'class')
            Matrix.C.train <- caret::confusionMatrix( factor(y.ouput), factor(data$V1[-indices])) #error de train
            Y.test.hat <- predict(y_hat, data[indices,], type = 'class')
            res2 <- Y.test.hat
            Matrix.C.test <- caret::confusionMatrix( factor(res2), factor(data$V1[indices]))
            acc.train[i] <- Matrix.C.train$overall['Accuracy']
            acc.test[i] <-  Matrix.C.test$overall['Accuracy']
        }
        gc()
        return( list(train=acc.train, test =acc.test) )
    }
}
train.SVM <- train.SVM.init(cv = 10, data = data, x)
error.SVM <- train.SVM(.1)#determinado por busqueda exaustiva en [.1, 20]
error.SVM <- as.data.frame(error.SVM)
error.SVM$modelo <- 'SVM'
############################################  busqueda en grid del parametro C de SVM el mejor es .1
# train.SVM <- train.SVM.init(cv = 10, data = data, x)
# c.tuning <- seq(0.1, 20, length=15)
# searh.grid.SVM <- mclapply(FUN = train.SVM,  X= c.tuning, mc.cores = (detectCores()-1))
# searh.grid.SVM
# z <- as.data.frame(searh.grid.SVM)
# z2 <- as.data.frame(apply(z, 2, mean))
# z.vis <- matrix(0, ncol = 2, nrow = 15)
# impares <-  seq(1,dim(z2)[1], by = 2)
# z.vis[,1 ] <- z2$`apply(z, 2, mean)`[impares]
# z.vis[,2 ] <- z2$`apply(z, 2, mean)`[impares+1]
# library(ggplot2)
# z.vis <- as.data.frame(z.vis)
# colnames(z.vis) <- c('Acc.Train', 'Acc.Test')
# z.vis$Indice.Grid <- 1:length(c.tuning)
#saveRDS(z.vis, file ='tuning_SVM.rds')
z.vis <-readRDS(file='tuning_SVM.rds')
p2 <- ggplot(z.vis, aes(x=Acc.Train, y=Acc.Test , color = factor(Indice.Grid)))+
    geom_point() +  theme(legend.title = element_blank()) +theme_minimal() +
    ggtitle('Precisión promedio por valor del grid (10-fold): SVM')
p2 <- ggplotly(p2) #distro en bayes
p2
############# vamos por el arbol################################
train.arbol.init <- function(cv = 10, data, grid)
{
    library(rpart)
    set.seed(0)
    bloque <- round(dim(data)[1]/cv)
    breaks <- c(seq(1, dim(data)[1], by = bloque  ), dim(data)[1])
    acc.train <- rep(0, cv)
    acc.test <- rep(0, cv)
    grid <- grid
    function(x){
        for (i in 1:(length(breaks)-1))
        {
            indices <- breaks[i]:breaks[i+1]
            test <- data[ indices,   ]
            train <- data[ -indices,   ]
            y_hat <- rpart(factor(V1)  ~ ., data = data[-indices,], method='class', 
                           control = rpart.control(cp = .01, minsplit= grid[x, 'minsplit' ],
                                                   maxdepth=grid[x, 'maxdepth' ]))
            y.ouput <-  predict(y_hat, data[-indices, ], type = 'class')
            Matrix.C.train <- caret::confusionMatrix( factor(y.ouput), factor(data$V1[-indices])) #error de train
            Y.test.hat <- predict(y_hat, data[indices,], type = 'class')
            res2 <- Y.test.hat
            Matrix.C.test <- caret::confusionMatrix( factor(res2), factor(data$V1[indices]))
            acc.train[i] <- Matrix.C.train$overall['Accuracy']
            acc.test[i] <-  Matrix.C.test$overall['Accuracy']
        }
        return( list(train=acc.train, test =acc.test) )
    }
}
#tuning el arbol 
grid <- data.frame(maxdepth = rep(1:15, each= 4), minsplit = rep((1:4)*5, 15))
train.arbol <- train.arbol.init(cv = 10, data = data, grid = grid)
error.arbol <- train.arbol(60) #resultado del tuning
error.arbol <- as.data.frame(error.arbol)
error.arbol$modelo <- 'arbol'
#tuning el arbol 
# grid <- data.frame(maxdepth = rep(1:15, each= 4), minsplit = rep((1:4)*5, 15))
# train.arbol <- train.arbol.init(cv = 10, data = data, grid = grid)
# searh.grid.arbol <- mclapply(FUN = train.arbol,  X= 1:dim(grid)[1], mc.cores = (detectCores()-1))
# z <- as.data.frame(searh.grid.arbol)
# z2 <- as.data.frame(apply(z, 2, mean))
# z.vis <- matrix(0, ncol = 2, nrow = dim(grid)[1])
# impares <-  seq(1,dim(z2)[1], by = 2)
# z.vis[,1 ] <- z2$`apply(z, 2, mean)`[impares]
# z.vis[,2 ] <- z2$`apply(z, 2, mean)`[impares+1]
# library(ggplot2)
# z.vis <- as.data.frame(z.vis)
# colnames(z.vis) <- c('Acc.Train', 'Acc.Test')
# z.vis$Indice.Grid <- 1:dim(grid)[1]
#saveRDS(z.vis, file ='tuning_arbol.rds')
z.vis <- readRDS(file='tuning_arbol.rds')
p2 <- ggplot(z.vis, aes(x=Acc.Train, y=Acc.Test , color = factor(Indice.Grid)))+
    geom_point() +  theme(legend.title = element_blank()) +theme_minimal() +
    ggtitle('Precisión promedio por valor del grid (10-fold): Arboles ')
p2 <- ggplotly(p2) #distro en bayes
p2


###################### NO  reutilizamos nuestra implementación de adaboost#################
set.seed(0)
train.ada.init <- function(cv = 10, data, grid)
{
    set.seed(0)
    library(adabag)
    bloque <- round(dim(data)[1]/cv)
    breaks <- c(seq(1, dim(data)[1], by = bloque  ), dim(data)[1])
    acc.train <- rep(0, cv)
    acc.test <- rep(0, cv)
    grid <- grid
    data$V1 <- factor(data$V1)
    function(x){
        indices <- sample(1:dim(data)[1], round(dim(data)[1]*.9))
        y_hat <- boosting( V1~ ., data= data[-indices,], boos=FALSE, 
                           coeflearn = 'Zhu', mfinal = grid[x, 'mfinal'], 
                           control=rpart.control(
                               maxdepth= grid[x, 'maxdepth'] , 
                               minsplit = grid[x, 'minsplit' ], cp =0.01 ) ) 
        
        y.ouput <-  predict(y_hat, data[-indices, ], type = 'class')
        Matrix.C.train <- caret::confusionMatrix( factor(y.ouput$class), factor(data$V1[-indices])) #error de train
        Y.test.hat <- predict(y_hat, data[indices,], type = 'class')
        res2 <- Y.test.hat
        Matrix.C.test <- caret::confusionMatrix( factor(res2$class), factor(data$V1[indices]))
        acc.train[x] <- Matrix.C.train$overall['Accuracy']
        acc.test[x] <-  Matrix.C.test$overall['Accuracy']
        gc()
        return( list(train=acc.train, test =acc.test) )
    }
}
grid.ada <- data.frame(minsplit = rep(seq(5, 20, by=2),  each = 9 ),
                       maxdepth = rep(seq(5,  15, by =5 ), each =  3), 
                       mfinal = rep(c(50, 100, 500), 4) )
grid.ada <- unique(grid.ada)
train.ada <- train.ada.init(cv = 10, data = data, grid = grid.ada)
error.ada <- train.ada(24) #determinado con busqueda en grid
error.ada <- as.data.frame(error.ada)
error.ada$modelo <- 'ada'
##################tuning de ada boost ######################################
# searh.grid.ada <- mclapply(FUN = train.ada,  X= 1:dim(grid.ada)[1], mc.cores = (detectCores()-1))
# o.2 <- searh.grid.ada
# o.3 <-o.2
# error.ada <- unlist(searh.grid.ada)
# #error.ada <- as.data.frame(error.ada)
# z <- matrix(-1, ncol = 2, nrow = 72)
# for(i in 1:length(o.2))
# {
#     foo1 <- o.2[[i]]$train
#     foo1[foo1==0] <- NA 
#     z[i, 1] <- mean(foo1, na.rm = TRUE)
#     foo1 <- o.2[[i]]$test
#     foo1[foo1==0] <- NA 
#     z[i, 2] <- mean(foo1, na.rm = TRUE)
# }
# z <- as.data.frame(z)
# names(z) <- c('Acc.Train', 'Acc.Test')
# z$Indice.Grid <- 1:72
library(ggplot2)
z.vis <- z
#saveRDS(z.vis, file ='tuning_ada.rds')
readRDS(file='tuning_ada.rds')
p2 <- ggplot(z.vis, aes(x=Acc.Train, y=Acc.Test , color = factor(Indice.Grid)))+
    geom_point() +  theme(legend.title = element_blank()) +theme_minimal() +
    ggtitle('Precisión promedio por valor del grid (10-fold): Arboles ')
p2 <- ggplotly(p2) #distro en bayes
p2
###########3grafico para comparar los desempellos en el dataset de numeros ##########
gc()
library(ggpubr)
errores <- rbind(error.lm, error.lda, error.qda, error.nn, error.SVM, error.arbol, error.ada)
p1 <-ggplot(errores, aes(x = modelo, y = train, color = modelo)) + geom_boxplot() +geom_jitter(width = 0.2, alpha = 0.4, color ='navy')+
    theme_minimal() +xlab('') + ylab('Presición en conjunto de entrenamiento')+ 
    ggtitle('Desempeño comparado entre modelos en 10-Folds (dataset números)')
p2 <- ggplot(errores, aes(x = modelo, y = test, color = modelo)) + geom_boxplot() +geom_jitter(width = 0.2, alpha = 0.5, color = 'navy')+
    theme_minimal() +xlab('') + ylab('Presición en conjunto de prueba')+ 
    ggtitle('Desempeño comparado entre modelos en 10-Folds (dataset números)')
ggarrange(p1,p2)
