####################################
#####    J. Antonio Garcia #########
#####   jose.ramirez@cimat.mx ######
####################################
library(plotly)
library(MASS)
library(shiny)
library(knitr)
library(rmarkdown) 
#setwd('/home/fou/Desktop/MCE_CIMAT/Second/CienciaDeDatos/DWD')
#########################################
# Construccion del backend              # 
#########################################
#lectura de datos simulados previamente
n <- 20
load('datos.RData')
####################################
server <- function(input, output) {
  #############################markdown de SVM y DWD
  # output$SVM <- renderUI({  
  #   rmarkdown::render(input = "SVM.Rmd",
  #                     output_format = html_document(self_contained = TRUE),
  #                     output_file = 'SVM.html')  
  #   shiny::includeHTML('SVM.html') 
  # })
  
  #seleccion de datos simulados
  a <- reactive({
    d <- input$d #seleccion de datos
    I <- diag(rep(1, d))
    pos <- pos[, 1:d]
    neg <- neg[, 1:d]
    stack <- rbind(pos,neg)
    pos.mean <- apply(pos, 2, mean)
    neg.mean <- apply(neg, 2, mean)
    w <- ginv(I) %*% (pos.mean - neg.mean)
    w <- w/sum(w**2)**.5 #normalizamos el vector MDP
    X <- ginv(cov(pos)) %*% (pos.mean - neg.mean)
    X <- X/sum(X**2)**.5 #normalizamos el vector que define al frontera de Bayes
    #acos(sum(X*w))*360/(2*pi)
    stack$label <- 1
    stack$label[(n+1):(2*n)] <- -1 
    Y <- X
    Y[-c(1,2)] <- 0
    Y[1] <- X[2]
    Y[2] <- -X[1] #elegimos un vector perpenticular a X
    Y <- Y/sum(Y**2)**.5
    sum(Y*w)
    M <- cbind(X, Y,  w)
    pos.proyec <- as.matrix(pos)%*%M
    neg.proyec <- as.matrix(neg)%*%M
    pos.proyec <- as.data.frame(pos.proyec)
    neg.proyec <- as.data.frame(neg.proyec)
    pos.proyec$label <- '+1'
    neg.proyec$label <- '-1'
    proyec <- rbind(pos.proyec, neg.proyec)
    proyec$w <- proyec$V2 * (w[1]/w[2])
    b <- list(proyec, w )
    return(b)
  })
  
  
  #construccion de scaterplot del tab DWD
  output$puntos <- renderPlotly({
    datos <- a()
    proyec <- datos[[1]]
    w <- datos[[2]]
    fit <- lm(w ~ V2, data = proyec)
    pal <- c("purple", "orange")
    p1 <- ggplot(data = proyec, aes(x=V2, y=V3, color=label))+geom_point() +
      stat_function(fun = function(z){z*(w[1]/w[2])}, aes(colour = I('MDP')), size=1.5)+
      geom_hline(yintercept=0, aes(colour=I('red')),     show.legend = NA)+
      ggtitle('Proyección en la dirección opima de Bayes y MaxDataPiling') +
      theme_minimal()+  xlab('Bayes') + ylab('') + 
      scale_colour_manual(
        labels = c('-1', '+1', 'MDP'),
        values = c("purple", "orange", "green4")
      )+ theme(legend.title = element_blank())
    p2 <- ggplotly(p1) #distro en bayes
    p2 #puntos pegados a la frontera
  })
  
  
  
  #construccion de densidades del tab DWD
  output$Bayes <- renderPlotly({
    a <- a()
    proyec <- a[[1]]
    p2 <- ggplot(data = proyec, aes(x=V3, fill=label, colour=label))+geom_density()+
      geom_rug(sides="b")+ggtitle('Distribución en la dirección de Bayes') + theme_minimal()+
      xlab('Bayes') + ylab('') + 
      scale_fill_manual(  labels = c('-1', '+1'), values = c("purple", "orange"))+
      scale_color_manual(  labels = c('-1', '+1'), values = c("purple", "orange"))+
      theme(legend.title = element_blank())
    p2 <- ggplotly(p2) #distro en bayes
    p2
  })
  
  #construccion de proyecciones del tab DWD
  output$DWD <- renderPlotly({
    a <- a()
    proyec <- a[[1]]
    x <- subset(proyec, label == '+1')
    y <- subset(proyec, label == '-1')
    y$V3 <- abs(y$V3)#distancia al plano separador
    p <- plot_ly(y = ~x$V3, type = "box",
                 line = list(color = 'rgb(255,165,0)'),
                 name = "+1") %>%
      add_trace(y = ~y$V3,
                line = list(color = 'rgb(160,32,240)'),
                name = "-1"
                ) %>%
      layout(title = "Distancia a hiperplano MDP",
             xaxis = list(title = ' '),
                          yaxis = list(title = ''))
    
  })
}
