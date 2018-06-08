####################################
#####    J. Antonio Garcia #########
#####   jose.ramirez@cimat.mx ######
####################################
library(readxl)
library(mxmaps)
library(plotly)
library(ggplot2)
library(mxmaps)
library(wordcloud)
library(wordcloud2)
data("df_mxstate")
#####################################3conocimiento apriori
actividades <- c('Agricultura, cría y explotación de animales...',
  'Minería',
  'Generación, transmisión y distribución de energía eléctrica...',
  'Construcción', 'Industrias manufactureras','Comercio al por mayor','Comercio al por menor',
  'Transportes, correos y almacenamiento','Información en medios masivos',
  'Servicios financieros y de seguros', 'Servicios inmobiliarios y de alquiler de bienes',
  'Servicios profesionales, científicos y técnicos','Corporativos',
  'Servicios de apoyo a los negocios y manejo de desechos...',
  'Servicios educativos', 'Servicios de salud y de asistencia social',
  'Servicios de esparcimiento culturales y deportivos...',
  'Servicios de alojamiento temporal...',
  'Otros servicios excepto actividades gubernamentales',
  'Actividades legislativas, gubernamentales, ...','Otro')
setwd('/home/fou/Desktop/MCE/MCE/Second/EstadisticaMultivariada/ProyectoFinal/CONOCER_Data_Docs/CONOCER_Data_Docs')
#######3lectura de datos
encuesta <- read_excel(path = 'CONOCER.xlsx', sheet = 'CONOCER FINAL' )
encuesta$cuo01 <- as.character(encuesta$cuo01)
ocupaciones <- read.csv('A.csv')
ocupaciones <- unique(ocupaciones) #por si pegue dos veces el mismo archivo 
ocupaciones$Clave <- as.character(ocupaciones$Clave)
ocupaciones <- ocupaciones[ order(ocupaciones$Clave), ] #el rollo es que la pagina da diferentes ocupaciones para la misma clave
###########################  limpieza del conjunto de ocupaciones con descripcion
ocupaciones.clave <- unique(ocupaciones$Clave)
index <- rep(0, length(ocupaciones.clave))
############################  
for(i in ocupaciones.clave)
{
    indice <- which( i == ocupaciones$Clave  )[1]
    index[i] <- indice
}
ocupaciones <- ocupaciones[index, ]
encuesta <- merge(encuesta, ocupaciones, all.x =TRUE, by.x ='cuo01', by.y ='Clave' ) #agregamos la descripcion del empleo
encuesta$act.eco <- actividades[encuesta$p002]
################3renombramos los estados 
encuesta$estado <- factor(encuesta$estado)
indices <- c(1, 2, 3, 6, 7, 8, 6, 9, 10,15, 11, 13, 14, 16,17, 19,20, 21, 22, 23, 24, 26, 27, 28,
  29, 30, 31)# identificacion manual
levels(encuesta$estado) <-  df_mxstate$state_name[indices] 
encuesta$estado <- as.character(encuesta$estado)
indices <- match(encuesta$estado, df_mxstate$state_name) 
encuesta$region <- df_mxstate$region[indices]
length(unique(encuesta$cuo01)) #443 grupos
length(unique(encuesta$empresa)) #14498 empresas grupos por errores de tipeo se repiten algunas 
#saveRDS(encuesta, 'encuesta.rds')
################################# #Distribucion de los empleados por estado
foo1 <- as.data.frame(table(encuesta$region)) 
foo1 <- foo1[order(foo1$Freq), ]
foo1$value <- foo1$Freq
foo1$region <- as.character(foo1$Var1)
a <- mxstate_choropleth(foo1, title = "Partición de los empleados en la encuesta ", 
                   legend = 'Número', num_colors = 1, zoom = foo1$region ) 
a
a <- ggplotly(a)
a
##################  #Distribucion de las empresas IMPORTANTES 
foo2 <- as.data.frame(table(encuesta$empresa)) 
foo2 <- foo2[order(foo2$Freq, decreasing = TRUE), ]  # cuantas empresas 14498
foo2$Var1 <- as.character(foo2$Var1)
foo2$Var1 <-factor(foo2$Var1 ,  levels = foo2$Var1[order(foo2$Freq)])
#length(unique(encuesta$empresa))
ggplot(head(foo2, 20), aes(x=Var1, y = Freq)) +geom_bar( stat="identity", fill = I('#6576B5')) + 
    coord_flip() +  theme_minimal()  + xlab('Nombre d ela empresa') + ylab('No. de encuestados')+
    ggtitle('Top 20 de empresas en la encuesta')
############################################## Participacion de empresas por estados
library(dplyr)
encuesta %>% group_by( empresa, region) %>% summarise(value = 1) -> foo
foo <- foo[order(foo$empresa, foo$region),]
foo %>% group_by(region) %>% summarise(value =sum(value)) -> foo
foo$region
a <- mxstate_choropleth(foo, title = "Partición de empresas por estado en la encuesta ", 
                        legend = 'Número', num_colors = 1, zoom = foo1$region ) 
a
a <- ggplotly(a)
a #son casi iguales los mapas por que d ecada empresa se eligieron a lo mas 3 empleados
############################### word 
a <- as.data.frame(table(encuesta$p001))
a <- a[order(a$Freq, decreasing = TRUE),]
a$Var1 <-factor(a$Var1 ,  levels = a$Var1[order(a$Freq)])
p <- ggplot(head(a, 20), aes(x=Var1, y = Freq)) +geom_bar( stat="identity", fill = I('#6576B5')) + 
    coord_flip() +  theme_minimal()  + xlab('Actividad') + ylab('No. de encuestados')+
    ggtitle('Top 20 de actividades de las empresas')
p
sum(a$Freq<2) # actividades se registraron solo una vez
############################################# actividad economica
a <- as.data.frame(table(encuesta$p002))
a <- a[order(a$Freq, decreasing = TRUE),]
a$actividad <- actividades[a$Var1]
a$actividad <-factor(a$actividad ,  levels = a$actividad[order(a$Freq)])
p <- ggplot(head(a, 20), aes(x=actividad, y = Freq)) +geom_bar( stat="identity", fill = I('#6576B5')) + 
    coord_flip() +  theme_minimal()  + xlab('Actividad economica') + ylab('No. de encuestados')+
    ggtitle('Top 20 de actividades economicas de las empresas')
p
sum(a$Freq<2) # actividades se registraron solo una vez
#barplot(table(encuesta$p003)) #distribucion de los lugares de trabajo
############################### competencias futuras 
#En el último año ¿Cuáles fueron las  principales vacantes requeridas por esta empresa, institución u organización?
futuro <- encuesta[ , c('p00501','p00502','p00503','p00504','p00505','p00506','p00507') ]
futuro <- as.matrix(futuro)
i <- is.na(futuro)
futuro[i] <- NA
t <- as.data.frame(table(futuro))
t <- t[ order(t$Freq, decreasing = TRUE),]
names(t) <- c('word', 'freq')
wordcloud2(t[2:50,], color = "random-light", shuffle = FALSE )
############################### competencias a mediano plazo  
#¿Qué actividades o funciones productivas considera usted importantes en un mediano o largo plazo?
futuro <- encuesta[ , c('p01001', 'p01002', 'p01003') ]
futuro <- as.matrix(futuro)
i <- is.na(futuro)
futuro[i] <- NA
t <- as.data.frame(table(futuro))
t <- t[ order(t$Freq, decreasing = TRUE),]
names(t) <- c('word', 'freq')
wordcloud2(t[2:50,], color = "random-light", shuffle = FALSE )
################################################ riesgos laborales p020x
futuro <- encuesta[ , c('p01901','p01902','p01903','p01904','p01905','p01906') ]
futuro <- as.matrix(futuro)
i <- is.na(futuro)
futuro[i] <- NA
t <- as.data.frame(table(futuro))
t <- t[ order(t$Freq, decreasing = TRUE),]
names(t) <- c('word', 'freq')
wordcloud2(t[2:50,], color = "random-light", shuffle = FALSE )
############ nivel de estudios ideal para el puesto 
futuro <- encuesta[ , c('p020') ]
futuro <- as.matrix(futuro)
n <- dim(encuesta)[1]
t <- as.data.frame(table(futuro)/n)
t$futuro <- c('Ninguno', 'Preescolar', 'Primaria','Secundaria','Carrera técnica',
              'Bachillerato tecnológico', 'Profesional técnico (técnico medio)', 
              'Bachillerato general (preparatoria)', 'Técnico superior universitario', 
              'Normal', 'Profesional superior','Maestría' ,'Doctorado')
t$futuro <-factor(t$futuro ,  levels = t$futuro[order(t$Freq)])
names(t) <- c('word', 'freq')
p <- ggplot(t, aes(x=word, y = freq)) +geom_bar( stat="identity", fill = I('#6576B5')) + 
    coord_flip() +  theme_minimal()  + xlab('') + ylab('% de No. de encuestados')+
    ggtitle('Nivel educativo requerido para el puesto actual')
p
########################### Experiencia necesaria
futuro <- encuesta[ , c('p021') ]
futuro <- as.matrix(futuro)
n <- dim(encuesta)[1]
t <- as.data.frame(table(futuro)/n)
t$futuro <- c('Ninguna', 'Hasta un máximo de 1 mes', 'Más de 1 mes hasta un máximo de 3 meses',
              'Más de 3 meses hasta un máximo de 6 meses','Más de 6 meses hasta un máximo de 1 año',
              'Más de 1 año hasta un máximo de 2 años','Más de 2 años hasta un máximo de 4 años',
              'Más de 4 años hasta un máximo de 6 años','Más de 6 años hasta un máximo de 8 años',
              'Más de 8 años hasta un máximo de 10 años', 'Más de 10 años')
t$futuro <- factor(t$futuro,  levels = t$futuro[order(t$Freq)])
names(t) <- c('word', 'freq')
p <- ggplot(t, aes(x=word, y = freq)) +geom_bar( stat="identity", fill = I('#6576B5')) + 
    coord_flip() +  theme_minimal()  + xlab('') + ylab('% de No. de encuestados')+
    ggtitle('Experiencia requerida para el puesto actual')
p
############################# 11-18
names(encuesta)
variables <- 38:139
data <- encuesta[, variables ]
saveRDS(data, 'unfold.rds')
etiquetas <- encuesta$Descripcion
saveRDS(etiquetas, 'etiquetas.rds')
#sum(is.na(data)) # hay nulos 
#remove(list = ls())
#gc()
setwd('/home/fou/Desktop/MCE/MCE/Second/EstadisticaMultivariada/ProyectoFinal/CONOCER_Data_Docs/CONOCER_Data_Docs')
library(smacof)
data <- readRDS('unfold.rds')
etiquetas <- readRDS('etiquetas.rds')
colnames(data)
rownames(data)
maxi <- max(data)
data <- (maxi+1) - data
sum(data<0)
###############Representacion solicitada en 2 y 35 dimensiones
encuesta <- readRDS(file='encuesta.rds')
etiquetas <- readRDS('etiquetas.rds')
t1 <- Sys.time()
t1
p <- unfolding(data, ndim = 35, type = 'mspline' ) #2: 25 minutos, 10: 13mins, 20: 5 mins
#p <- unfolding(data, ndim = 35) #35: 3.5 minutos 
t1 <- Sys.time() -t1
print(t1)
gc()
t1
save(p,file= 'p2.rdata')
save(p, file ='/srv/shiny-server/sample-apps/CONOCER/p2.rdata')
##############################################
setwd('/srv/shiny-server/sample-apps/CONOCER/')
indices <- match( as.character(etiquetas), levels(etiquetas))
colores.grupo <- colors()[indices%% length(colors()) ]
#plot(p , label.conf.rows = list(label = FALSE), main ='MDS unfolding usando 2 dimensiones stress(0.1362)')
plot(p$conf.col , main ='MDS unfolding usando 35 dimensiones (stress 0.01678)',
     xlim=c(-.65,.9), pch=20, col = 'purple' , ylim = c(-.6,.6)  , cex =1. )
points(p$conf.row, col = nuevos.colores, cex=.5 , pch = 20)
text( p$conf.col  , labels = colnames(data), cex = .9, pos = sample(1:4, 102, replace = TRUE),
      col ='purple')
muestrita <- encuesta$Descripcion[colores.grupo=='black']
############################################### simulacion para numero de dimension en MDS
setwd("/home/fou/Desktop/MCE/MCE/Second/EstadisticaMultivariada/ProyectoFinal/CONOCER_Data_Docs/CONOCER_Data_Docs")
data <- readRDS('simulacion.rds')
set.seed(0)
dimensiones <- c(2,5,15,30,50) #aprovechando que dispongo de 7 hilos
library(parallel)
library(smacof)
t1 <- Sys.time()
stress <- mclapply(  X=dimensiones, FUN = function(i){
            sub <- data[sample(1:17250, 8625), ] # la mitad del tamaño 
            modelo <- unfolding(sub, ndim = i)
            gc()
            return(modelo$stress)
    }, mc.cores = 2)
t1 <- Sys.time() -t1
t1 # 1.891209 mins
numeros <- unlist(stress)
numeros <- as.data.frame(t(matrix(numeros, byrow = TRUE, ncol = 5)))
numeros$dimension <- dimensiones
numeros$simulacion <- '1'
gc()
##################################################
t1 <- Sys.time()
stress2 <- mclapply(  X=dimensiones, FUN = function(i){
    sub <- data[sample(1:17250, 8625), ] # la mitad del tamaño 
    modelo <- unfolding(sub, ndim = i)
    gc()
    return(modelo$stress)
}, mc.cores = 2)
t1 <- Sys.time() -t1
t1 # 32 sec
numeros2 <- unlist(stress2)
numeros2 <- as.data.frame(t(matrix(numeros2, byrow = TRUE, ncol = 5)))
numeros2$dimension <- dimensiones
numeros2$simulacion <- '2'
simula <- rbind(numeros, numeros2)
library(ggplot2)
ggplot(simula, aes(x= dimension, y = V1, color=simulacion ))+ geom_line()+ theme_minimal()+
    ylab('Stress')+ggtitle('Simulación para elección de número de dimensión')+xlab('no. dimensiones')
ls() ##borramos todo menos los datos
##################simulacion de PCA 
m <- 102
n <- 17250
library(parallel)
library(abind)
set.seed(0)
out <- mclapply(X=1:10, FUN = function(i){
    pc <- matrix( sample(1:5, m*n, replace =TRUE), ncol = m)
    descom <- eigen(cor(pc))
    descom <- matrix(descom$values, ncol = 1) 
    descom <- as.data.frame(descom)
    descom$simulacion <-(i)
    descom$n <- 1:102
    gc()
    return(descom)
})
out <- abind(out, along = 1)
out <- as.data.frame(out)
head(out)
apply(out, 2, class)
out$V1
library(ggplot2)
empiricos <- data
#head(empiricos)
descom <- eigen(cor(empiricos))
descom <- matrix(descom$values, ncol = 1) 
descom <- as.data.frame(descom)
descom$simulacion <- 'Muestral'
descom$n <- 1:102
head(descom)
#out$simulacion <- factor(out$simulacion)
head(out)
p1 <- ggplot(out, aes(x = n, y = V1, color = factor(simulacion)  )) + geom_point() +
    theme_minimal() + ggtitle('Valores propios de 100 conjuntos de datos simulados') + 
    xlab('n-esimo valor propio')+ylab('Valor')
p1
out <- rbind(out, descom)
p2 <- ggplot(out, aes(x = n, y = V1, color = factor(simulacion)  )) + geom_point() +
    theme_minimal() + ggtitle('Valores propios de los datos simulados y los muestrales ') + 
    xlab('n-esimo valor propio')+ylab('Valor')
p2
p3 <- ggplot(out, aes(x = n, y = V1, color = factor(simulacion)  )) + geom_point() +
    theme_minimal() + xlim(c(3,25))  + ylim(c(0.5,2.5))+ xlab('n-esimo valor propio')+ylab('Valor')+
    ggtitle('Valores propios de los datos simulados y los muestrales (acercamiento) ')
p3
library(ggpubr)
ggarrange(p1, p2, p3, ncol = 2, nrow = 2)
#############################con componentes 15 basta , generacion de indice 
descom <- prcomp(cor(empiricos))
summary(descom)
a <- descom$sdev
L <- descom$rotation[,1:2]
(round(L,1))
gente <- as.matrix(empiricos) %*% L
grupos <- encuesta$cuo01
indices <- match(grupos, unique(grupos))
nuevos.colores <- colors()[indices]
df <- as.data.frame(gente)
df$color <- nuevos.colores
ggplot(df, aes(x=PC1, y=PC2, color=color, alpha = 0.00001))+geom_point(size = .5) +
    theme_minimal()+ theme(legend.position="none") + xlab('Primer componente principal')+
    ylab('Segunda componente principal')+ggtitle('Resultado obtenido usando PCA')
