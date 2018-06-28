setwd('C:\\Users\\fou-f\\Desktop\\MCE\\Second\\CienciaDeDatos\\tarea2')
genes <- read.csv('gene_expression_2classes.csv', stringsAsFactors = FALSE, header=FALSE)
pacientes <- t(as.matrix(genes))
pacientes <- as.data.frame(pacientes)
S <- cor(pacientes)
e.d <- eigen(S) #realizamos PCA sobre la matriz de correlacion
cumsum(e.d$values)[1:5]/sum(e.d$values)#checamos la varianza explicada por las primeras 5componentes
datos.rotados <- (as.matrix(scale(pacientes)))%*%e.d$vectors #rotamos los datos iniciales escalados
set.seed(0)
grupos <- kmeans(x = (pacientes), centers = 2, algorithm = 'Lloyd', nstart = 10 )
plot(datos.rotados[,1], datos.rotados[,2], col= grupos$cluster, pch =20)
datos.rotados <- as.data.frame(datos.rotados)
clasificacion <- list('sano'=2, 'enfermo'=1)
clase <- unlist(names(clasificacion[grupos$cluster]))
library(ggplot2) #graficamos la proyeccion en las dos primeras componentes
ggplot(datos.rotados, aes(x=V1, y = V2,
                          color=clase)) +
  geom_point()+
  ggtitle('Kmeans sobre PCA (datos centrados)') + theme_minimal() + 
  xlab('Primer componente principal') +
  ylab('Segunda componente principal') 
cumsum(e.d$values[1:2])/sum(e.d$values)# pese a que la primer componente solo capta 8% de la variacion
distancia <- scale((as.matrix(pacientes))) #buscamos cluster usando metodos jerarquicos
arbol <- hclust(dist(distancia, method = 'euclidean' ), method ="ward.D")
plot(arbol)
labels <- cutree(arbol, k=2)
library(dendextend) #pintamos el arbol con los clusters
arbol.feliz <- arbol %>% as.dendrogram %>%
  set("branches_k_color", k=2, c('purple','green4')) %>% set("branches_lwd", 1.2) %>%
  set("labels_colors", c(rep('purple',20), rep('green4',20))) %>%
  set("labels_cex", c(.9)) 
plot(arbol.feliz)
################inciso 2
R <- cor((as.matrix(pacientes)))
a <- heatmap(t(as.matrix(genes))) #el heatmap parece reconocer grupos de genes parecidos 
image(t(as.matrix(genes))[a$rowInd, a$colInd]) #extraemos del heatmap los que parecen forman el grupo significativo 
columnas <- tail(a$colInd, 111)
image(R)
#realizamos cluster jerarquico con matriz de similaridades inducida por la matriz
# de correlacion
R_dist <- diag(1000)-R
distancias <- as.dist(R_dist) 
h <- hclust(distancias, method = 'ward.D')
plot(h)
labels<- cutree(h, h = 10)
freq <- as.data.frame(table(labels))
#decoramos el arbol
otro.arbol.feliz <- h %>% as.dendrogram %>%
  set("branches_k_color", k=2, c('purple','green4')) %>% set("branches_lwd", 1.2) %>%
  set("labels_colors", c(rep('purple',freq$Freq[2]), 
                         rep('green4',freq$Freq[1]))) %>%
  set("labels_cex", c(.1)) 
plot(otro.arbol.feliz) 
genes <- labels[labels==2] #extraemos lel segundo grupo del heatmap
columnas
genes.columnas <- gsub("V", "", names(genes))#extraemos solo el numero de columna
genes.columnas <- as.numeric(genes.columnas) 
comparacion <- as.matrix(table(sort(genes.columnas), sort(columnas)))#lo comparo contra el grupo de menor cardinalidad encontrado
image(diag(111)-comparacion) ##color uniforme, son identicos los conjuntos
table(labels)
genes.columnas
######### Analisis analogo sin los genes importantes
test <- pacientes[,-(genes.columnas)]
S <- cor(test)
e.d <- eigen(S)
cumsum(e.d$values)[1:5]/sum(e.d$values)
datos.rotados <- (as.matrix(scale(test)))%*%e.d$vectors 
grupos <- kmeans(x = scale(test), centers = 2, algorithm = 'Lloyd' )
plot(datos.rotados[,1], datos.rotados[,2], col= grupos$cluster, pch =20)
plot(datos.rotados[,1], datos.rotados[,2], 
     col= c(rep('red',20), rep('blue',20)), pch =20)
datos.rotados <- as.data.frame(datos.rotados)
clasificacion <- list('sano'=2, 'enfermo'=1)
clase <- unlist(names(clasificacion[c(rep(2,20), rep(1,20))]))
library(ggplot2)
ggplot(datos.rotados, aes(x=V1, y = V2,
                          color=clase)) +
  geom_point()+
  ggtitle('Kmeans sobre PCA (datos centrados) sin genes importantes') + theme_minimal() + 
  xlab('Primer componente principal') +
  ylab('Segunda componente principal') 
cumsum(e.d$values[1:2])/sum(e.d$values)# pese a que la primer componente solo capta 8% de la variacion
