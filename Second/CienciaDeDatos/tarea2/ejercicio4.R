directorio <-'C:\\Users\\fou-f\\Desktop\\MCE\\Second\\CienciaDeDatos\\tarea2\\movie_reviews\\txtfiles'
setwd(directorio)
dir()
#primero veamos las palabras de las opiniones negativas 
#para contrastarlas con las palabras de las opiniones positivas y lograr
#realzar el contraste entre ambos grupos quitando las palabras que aparecen en 
#ambas bolsas de palabras
preprocesamiento <- function(x)
{
  #esta funcion tiene como finalidad extraer todas las palabras
  #regresa un dataset con las frecuencias acumuladas por palabras en el corpus
  # y una matriz de treminos-documentos
  # x (path): path en donde se encuentran lojados en disco duro los documentos
  library(tm)
  negativos <- Corpus(DirSource(x,recursive=TRUE),
               readerControl=list(language="en_US"))
  #preprosamiento
  negativos <- tm_map(negativos,stripWhitespace)
  negativos <- tm_map(negativos,removeNumbers) #tal vez los numero sean inportantes
  negativos <- tm_map(negativos,content_transformer(tolower))
  negativos <- tm_map(negativos,removePunctuation) #talvez los emoticones tambien sean importantes
  negativos <- tm_map(negativos,removeWords,stopwords("english"))
  negativos <- tm_map(negativos,stemDocument)
  ## obtiene matriz de terminos
  matriz.neg <- TermDocumentMatrix(negativos)#,control=list(minDocFreq=100))
  m.neg <- as.matrix(matriz.neg)
  v <- sort(rowSums(m.neg),decreasing=TRUE)
  d.neg <- data.frame(word = names(v),freq=v)
  d.neg <- d.neg[order(d.neg$word),]
  return(list(frecuencias= d.neg, mtd=m.neg ))
}
#######################
negativos <- preprocesamiento(dir()[1]) #preprocesamos los documentos marcados como negativos 
positivos <- preprocesamiento(dir()[2]) #preproseamiento de los documentos marcados como positivos
# identificamos las palabras comunes a ambos conjuntos  
comunes <- merge(negativos[['frecuencias']],positivos[['frecuencias']]  , by ='word')
#checamos la correlacion de las frecuencias de las palabras en ambos conjuntos 
#primero visualmente
library(ggplot2)
ggplot(comunes, aes(x = freq.x, y = freq.y, alpha=.0001 ))+
  geom_point()+theme_minimal()+ theme(legend.position="none")+
  xlab('Conjunto de palabras en la bolsa "negativas"')+
  ylab('Conjunto de palabras en la bolsa "positivas"')+
 ggtitle('Frecuencia de palabas comunes en ambos conjuntos de reseñas ')
#parece indicar correlación positiva
#hacemos un test de significancia 
cor.test(comunes$freq.x, comunes$freq.y, method = 'pearson', alternative='two.sided')
#como el p-value es casi cero rechazamos la hipotesis nula de no correlacion 
#asumimos quue la frecuencia de las palabras es aproximadamente la misma 
#por lo que las descartamos del analisis pues solo incrementarian la dimención de la tarea y no 
#ayudan a diferenciar el sentimiento
#################################################
todos <- preprocesamiento(dir()) #todos los datos 
utiles <- todos[['mtd']]
#primer filtro quitamos las palabras que no aportan información 
utiles <- utiles[!(row.names(utiles) %in% as.character(comunes$word)),   ]
df.m <- as.data.frame(t(utiles))
d <- todos[['frecuencias']]
d <- d[!(d$word %in% as.character(comunes$word)),]
#procedemos a eliminar las palabras con menor frecuencia
poquitas <- d$word[d$freq>5] #determinamos las palabras comunes y con una frecuencia mayor a 6
#library(wordcloud)
wordcloud(comunes$word, comunes$freq.x)
utiles <- utiles[(row.names(utiles) %in% as.character(poquitas)),   ]
library(dplyr) #fiiltramode ugual manera la tabla de frecuuencias
df.m <- as.data.frame(t(utiles))
v <- sort(rowSums((utiles)),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
#nos quedamos con 443 palabras unicamente
##########seccion un usar PCA para distinguir las recomendaciones postivas de las negativas
#PCA sobre matriz de terminos 
R <- cor(df.m)
det(R)
pca <- eigen(R)
cumsum(pca$values)[1:20]/sum(pca$values) #las dos primeras componentes solo explican el 2.8% de la varianza total
rotacion <- pca$vectors
datos.rotados <- as.matrix(df.m)%*%rotacion #rotamos los datos
#graficamos las reseñas en el espacio de componentes principales
datos.rotados <- as.data.frame(datos.rotados)
datos.rotados$sentimiento <- 'negativo'
datos.rotados$sentimiento[501:1000] <- 'Positivo'
p1 <- ggplot(subset(datos.rotados, sentimiento=='negativo'), 
       aes(x = V1, y = V2, color=sentimiento, alpha=.01 ))+
  geom_point()+theme_minimal()+
  xlab('Primer componente principal')+
  ylab('Segundo componente principal')+
  ggtitle('Reseñas negativas proyectadas en las dos primeras componentes principales')
p2 <- ggplot(subset(datos.rotados, sentimiento=='Positivo'), 
       aes(x = V1, y = V2, color=sentimiento, alpha=.01 ))+
  geom_point()+theme_minimal()+
  xlab('Primer componente principal')+
  ylab('Segundo componente principal')+
  ggtitle('Reseñas positivas proyectadas en las dos primeras componentes principales') 
p3 <- ggplot(datos.rotados, aes(x = V1, y = V2, color=sentimiento, alpha=.0001 ))+
  geom_point()+ theme_minimal()+
  xlab('Primer componente principal')+
  ylab('Segundo componente principal')+
  ggtitle('Todas las reseñas proyectadas en las dos primeras componentes principales') +
  ylim(c(-1/10,1/10)) +xlim(c(-1/10,1/10))
library(ggpubr)
ggarrange(p1, p2, p3, ncol = 1, nrow = 3)
datos.rotados <- datos.rotados%>% 
  mutate(etiquetaNegativo = !( (V1>=0) & (V2<=0)))
table(datos.rotados$sentimiento, datos.rotados$etiquetaNegativo)
############segundo inciso 
# procedemos de manera analoga utilizando los datos que ya calculamos, 
#primero determinamos las palabras en ambos conjuntos pero agrupados de 5 en 5
#primero agrupamos las palabras en comun resultaron ser las mismas
# matriz.neg <- negativos[['mtd']]
# matriz.neg <- as.data.frame(t(matriz.neg))
# k <- 5
# matriz.neg$grupo <- factor(rep(1:(dim(matriz.neg)[1]/k), each=k))
# library(dplyr)
# test <- matriz.neg %>% group_by(grupo) %>% summarise_all(funs(sum))
# test <- as.matrix(as.numeric(test[,2:15622]))
# v <- apply(test[,-1], 2, sum)
# d.neg <- data.frame(word = names(v),freq=v)
# matriz.pos <- positivos[['mtd']]
# matriz.pos <- as.data.frame(t(matriz.pos))
# matriz.pos$grupo <- factor(rep(1:(dim(matriz.pos)[1]/k), each=k))
# test2 <- matriz.pos%>%group_by(grupo)%>%summarise_all(funs(sum))
# v <- apply(test2[,-1], 2, sum)
# d.pos <- data.frame(word = names(v),freq=v)
# #por fin obtenemos las palabra comunes
# comunes <- merge(d.pos, d.neg  , by ='word') #las mismas 9477 palabras
#################################################
todos <- preprocesamiento(dir()) #todos los datos 
utiles <- todos[['mtd']]
#primer filtro quitamos las palabras que no aportan información 
utiles <- t(utiles)
k <- 5 #numero de agrupacion 
utiles <- as.data.frame(utiles)
utiles$grupo <- factor(rep(1:(dim(utiles)[1]/k), each=k))
test <- utiles %>% group_by(grupo) %>% summarise_all(funs(sum))
dim(test)
head(test[,1:5])
v <- apply(test[,-1], 2, sum)
d <- data.frame(word = names(v),freq=v)
utiles.fino <- test[, !(colnames(test) %in% as.character(comunes$word))   ]
#procedemos a eliminar las palabras con menor frecuencia
poquitas <- d$word[d$freq>5] #determinamos las palabras comunes y con una frecuencia mayor a 6
utiles.fino <- utiles.fino[, (colnames(utiles.fino) %in% as.character(poquitas))   ]
df.m <- as.data.frame(utiles.fino)
v <- apply(utiles.fino[,-1], 2, sum)
d <- data.frame(word = names(v),freq=v)
#nos quedamos con 442 palabras unicamente
##########seccion un usar PCA para distinguir las recomendaciones postivas de las negativas
#PCA sobre matriz de terminos 
R <- cor(df.m)
det(R)
pca <- eigen(R)
cumsum(pca$values)[1:20]/sum(pca$values) #las dos primeras componentes solo explican el 4.7% de la varianza total
rotacion <- pca$vectors
datos.rotados <- as.matrix(df.m)%*%rotacion #rotamos los datos
#graficamos las reseñas en el espacio de componentes principales
datos.rotados <- as.data.frame(datos.rotados)
datos.rotados$sentimiento <- 'negativo'
datos.rotados$sentimiento[101:200] <- 'Positivo'
p1 <- ggplot(subset(datos.rotados, sentimiento=='negativo'), 
             aes(x = V1, y = V2, color=sentimiento, alpha=.01 ))+
  geom_point()+theme_minimal()+
  xlab('Primer componente principal')+
  ylab('Segundo componente principal')+
  ggtitle('Reseñas negativas proyectadas en las dos primeras componentes principales')
p2 <- ggplot(subset(datos.rotados, sentimiento=='Positivo'), 
             aes(x = V1, y = V2, color=sentimiento, alpha=.01 ))+
  geom_point()+theme_minimal()+
  xlab('Primer componente principal')+
  ylab('Segundo componente principal')+
  ggtitle('Reseñas positivas proyectadas en las dos primeras componentes principales') 
p3 <- ggplot(datos.rotados, aes(x = V1, y = V2, color=sentimiento, alpha=.0001 ))+
  geom_point()+ theme_minimal()+
  xlab('Primer componente principal')+
  ylab('Segundo componente principal')+
  ggtitle('Todas las reseñas proyectadas en las dos primeras componentes principales') 
ggarrange(p1, p2, p3, ncol = 1, nrow = 3)
datos.rotados <- datos.rotados%>% 
  mutate(etiquetaNegativo = !( (V1<=0) & (V2<=0)))
table(datos.rotados$sentimiento, datos.rotados$etiquetaNegativo)
