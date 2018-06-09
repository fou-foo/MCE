####################################
#####    J. Antonio Garcia #########
#####   jose.ramirez@cimat.mx ######
####################################
## Evaluacion del metodo con datos reales
############### cancer
setwd('/home/fou/Desktop/MCE/Second/CienciaDeDatos/DWD1')
dir()
cancer <- read.csv('cancer.csv', header = FALSE)
grid <- seq(0.3162278, 31, length =5  )
grid <- grid**2
cancer$V1 <- NULL
y <- cancer$V2
cancer$V2 <- NULL
data <- cancer
train.lm.cv <- function(x, data, y  )
{
    data <- data
    y <- y 
    cv <- 10
    #funcion para evaluar por cv lm
    set.seed(0)
    #data <- data[ sample(1:dim(data)[1], dim(data)[1]),]
    p <- dim(data)[2]
    bloque <- round(dim(data)[1]/10)
    breaks <- c(seq(1, dim(data)[1], by = bloque  ), dim(data)[1])
    error <- matrix( rep(1, 10*5 ), ncol = 5  )
    colnames(error) <- c('MDP', 'RLR', 'SVM', 'DWD', 'Adaboost')
    library(MASS)
    function(x){
        for (i in 1:(length(breaks)-1))
        {
            indices <- breaks[i]:breaks[i+1]
            test <- data[ indices,   ]
            train <- data[ -indices,   ]
            y.test <- y[indices]
            y.train <- y[-indices]
            ##############evaluacion de MDP 
            pos.mean <- apply(train[ y.train=='B',  ], 2, mean)
            neg.mean <- apply(train[ y.train=='M',  ], 2, mean)
            X <- ginv(cov(train)) %*% (pos.mean - neg.mean)
            MDP <- X/sum(X**2)**.5 #normalizamos el vector 
            test.predic <- as.matrix(test)%*%MDP
            y_hat <-  ifelse(test.predic>=0, 1, -1)
            error.MDP <- sum(diag(table( y_hat, y.test)))/dim(test)[1]
            error[i, 'MDP'] <- error.MDP
            gc()
            ####################evaluacion de RLR
            library(glmnet)
            RLR <- glmnet(x = as.matrix(train), y = y.train , family = "binomial")
            #plot(RLR)
            y.RLR <- predict(RLR, newx = as.matrix(test), type = "class", s = x) #el paper dice que usaron este valor de lambda
            e <- table( y.RLR, y.test)
            error.RLR <- 1-sum(diag(e))/dim(test)[1]
            error[i, 'RLR'] <- error.RLR
            gc()
            ########################3 evaluacion de svm
            library(e1071)
            svm <- svm( y.train ~ ., data = cbind(train, y.train), kernel ='linear', cost = x)
            y.hat <- predict(svm, test)
            error.svm <- 1 - sum(diag(table(y.hat, y.test)))/dim(test)[1]
            error[i, 'SVM'] <- error.svm
            gc()
            #################### evaluacion de DWD
            library(kerndwd)
            kern <- vanilladot()
            dwd <- kerndwd(x=train, y.train, kern, lambda = x)
            y.hat <- predict.kerndwd(dwd,  kern=kern, newx = as.matrix(test), x = as.matrix(train), type='class')
            error.DWD <-  1-sum(diag(table( y.test, y.hat)))/dim(test)[1]
            error[i, 'DWD'] <- error.DWD
            gc()
            library(fastAdaboost)
            ada <- adaboost( y.train ~ . , data=cbind(train, y.train),  round(x)+1)   
            y.hat <- predict(ada, newdata = test, type ='class') 
            error.ada <- 1 - sum(diag(table(y.hat$class, y.test)))/dim(test)[1]
            error[i, 'Adaboost'] <- error.ada
        }
        medias <- apply(error, 2, mean)
        sds <- apply(error, 2, sd)
        res <- rbind(medias, sds)
        res <- as.data.frame(t(res))
        res$parametro <- as.character(x)
        gc()
        return( res )
    }
}
train.cv.10 <- train.lm.cv(data=cancer, y=y)
library(parallel)
simula.cancer <- mclapply(X=grid, FUN=train.cv.10, mc.cores = 4)
library(abind)
a <- abind(simula.cancer, along = 1)
a2 <- data.frame(metodo = row.names(a), media = as.numeric(a[,1]), sd = as.numeric(a[,2]), parametro = a[,3])
ggplot(a2, aes(x=as.numeric(parametro), y=media, colour=metodo)) + 
    #geom_errorbar(aes(ymin=media-sd, ymax=media+sd), width=.1) +
    geom_line() + geom_point() + scale_x_continuous(trans='log10') + theme_minimal()+ xlab('parametro')+
    ylab('Error de test') + ggtitle('Desempeño en Wisconsin Diagnostic Breast Cancer')
ggplot(a2, aes(x=as.numeric(parametro), y=media, colour=metodo)) + 
    #geom_errorbar(aes(ymin=media-sd, ymax=media+sd), width=.1) +
    geom_line() + geom_point() + scale_x_continuous(trans='log10') + theme_minimal()+ xlab('parametro')+
    ylab('Error de test') + ggtitle('Desempeño en Wisconsin Diagnostic Breast Cancer') +ylim(c(0,.1))
########################################################
####################### datos genes
setwd('/home/fou/Desktop/MCE/Second/CienciaDeDatos/DWD1')
dir()
genes <- read.csv('/home/fou/Desktop/TCGA-PANCAN-HiSeq-801x20531/data.csv', header = TRUE)
genes$X <- NULL
y <- read.csv('/home/fou/Desktop/TCGA-PANCAN-HiSeq-801x20531/labels.csv', header = TRUE)
table(y$Class)
#################3 muestreo por columnas
index <-  y$Class %in% c('COAD', 'LUAD')
set.seed(123)
datos <- genes[index, sample(1:20531,456*4) ]
y <- y$Class[index]
table(y)
######################3 muestreo por reglones
index <-  sample(1:219, 136)
datos <- datos[index, ]
y <- y[index]
table(y)
write.csv(datos, 'subgenes.csv')
write.csv(y, 'sub_y.csv')
genes <- read.csv('subgenes.csv')
data <- genes
y <- read.csv('sub_y.csv')
y <- y$x
data <- genes
grid <- seq(.1, 3, length = 10  )
train.lm.cv <- function(x, data, y  )
{
    data <- data
    y <- y 
    cv <- 5
    #funcion para evaluar por cv lm
    set.seed(0)
    #data <- data[ sample(1:dim(data)[1], dim(data)[1]),]
    p <- dim(data)[2]
    bloque <- round(dim(data)[1]/cv)
    breaks <- c(seq(1, dim(data)[1], by = bloque  ), dim(data)[1])
    error <- matrix( rep(1, 10*5 ), ncol = 5  )
    colnames(error) <- c('MDP', 'RLR', 'SVM', 'DWD', 'Adaboost')
    library(MASS)
    function(x){
        for (i in 1:(length(breaks)-1))
        {
            gc()
            indices <- breaks[i]:breaks[i+1]
            test <- data[ indices,   ]
            train <- data[- indices,   ]
            y.test <- y[indices]
            y.train <- y[-indices]
            ##############evaluacion de MDP 
            pos.mean <- apply(train[ y.train=='COAD',  ], 2, mean)
            neg.mean <- apply(train[ y.train=='LUAD',  ], 2, mean)
            X <- ginv(cov(train)) %*% (pos.mean - neg.mean)
            MDP <- X/sum(X**2)**.5 #normalizamos el vector 
            test.predic <- as.matrix(test)%*%MDP
            y_hat <-  ifelse(test.predic>=0, 1, -1)
            error.MDP <- sum(diag(table( y_hat, y.test)))/dim(test)[1]
            error[i, 'MDP'] <- error.MDP
            gc()
            ####################evaluacion de RLR
            library(glmnet)
            RLR <- glmnet(x = as.matrix(train), y = y.train , family = "binomial")
            #plot(RLR)
            y.RLR <- predict(RLR, newx = as.matrix(test), type = "class", s = x) #el paper dice que usaron este valor de lambda
            e <- table( y.RLR, y.test)
            error.RLR <- 1-sum(diag(e))/dim(test)[1]
            error[i, 'RLR'] <- error.RLR
            gc()
            ########################3 evaluacion de svm
            library(e1071)
            svm <- svm( y.train ~ ., data = cbind(train, y.train), kernel ='linear', cost = x)
            y.hat <- predict(svm, test)
            error.svm <- 1 - sum(diag(table(y.hat, y.test)))/dim(test)[1]
            error[i, 'SVM'] <- error.svm
            gc()
            #################### evaluacion de DWD
            library(kerndwd)
            kern <- vanilladot()
            dwd <- kerndwd(x=train, y.train, kern, lambda = x)
            y.hat <- predict.kerndwd(dwd,  kern=kern, newx = as.matrix(test), x = as.matrix(train), type='class')
            error.DWD <-  1-sum(diag(table( y.test, y.hat)))/dim(test)[1]
            error[i, 'DWD'] <- error.DWD
            gc()
            library(fastAdaboost)
            ada <- adaboost( y.train ~ . , data=cbind(train, y.train),  round(x)+1)   
            y.hat <- predict(ada, newdata = test, type ='class') 
            error.ada <- 1 - sum(diag(table(y.hat$class, y.test)))/dim(test)[1]
            error[i, 'Adaboost'] <- error.ada
        }
        medias <- apply(error, 2, mean)
        sds <- apply(error, 2, sd)
        res <- rbind(medias, sds)
        res <- as.data.frame(t(res))
        res$parametro <- as.character(x)
        gc()
        return( res )
    }
}
train.cv.10 <- train.lm.cv(data=genes, y=y)
library(parallel)
simula.cancer <- lapply(X=grid, FUN=train.cv.10)
library(abind)
a <- abind(simula.cancer, along = 1)
a2 <- data.frame(metodo = row.names(a), media = as.numeric(a[,1]), sd = as.numeric(a[,2]), parametro = a[,3])
library(ggplot2)
ggplot(a2, aes(x=as.numeric(as.character(parametro)), y=media, colour=metodo)) + 
    #geom_errorbar(aes(ymin=media-sd, ymax=media+sd), width=.1) +
    geom_line() + geom_point()  + theme_minimal()+ xlab('parametro')+
    ylab('Error de test') + ggtitle('Desempeño en genes')
ggplot(a2, aes(x=as.numeric(parametro), y=media, colour=metodo)) + 
    geom_errorbar(aes(ymin=media-sd, ymax=media+sd), width=.1) +
    geom_line() + geom_point()  + theme_minimal()+ xlab('parametro')+
    ylab('Error de test') + ggtitle('Desempeño en genes')
