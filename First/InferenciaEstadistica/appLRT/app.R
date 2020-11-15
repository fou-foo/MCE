library(shiny)
library(ggplot2)  # for the diamonds dataset
library(plotly)
library(lmtest)
library(shinydashboard)
library(reshape2)
library(MASS)
library("microbenchmark")
library(shinythemes)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("armadillo.cpp")
sourceCpp("armadillo_noparalelo.cpp")
## J.Antonio Garcia jose.ramirez@cimat.mx
###################
# xdim <- 2
# ydim <- 3
# 
# # Number of Time Steps
# tdim <- 2
# 
# # Generation of 3D Arrays
# 
# # Initial Neighborhoods
# a <- array(0:1, dim = c(xdim, ydim, tdim))
# 
# # Result storing matrix
# res <- array(0:1, dim = c(xdim, ydim, tdim))
# 
# 
# krit <- function(n){
#   mat <- matrix(1,n,1)
#   return(mat)
# }
# cppFunction('
#             // Added krit as a function pass
#             NumericMatrix tes1(int n, Rcpp::Function krit){
#             NumericMatrix test = krit(n+1);
#             return(test);
#             }
#             ')
n <- 1000
turnout1 <- data.frame(income <- rnorm(n),
                       educate <- rnorm(n, 2,3),
                       age <- rnorm(n, -1,20),
                       vote <- sample(0:1, n, replace = TRUE)
)
mY <- matrix(turnout1$vote)
mX <- cbind(1,
            turnout1$income,
            turnout1$educate,
            turnout1$age
)

Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

#######################
mtcars 
carros <- mtcars[,c("mpg", "disp","drat")]
modelo.actual<-lm(mpg ~  disp+drat, data = carros)
carros$z_hat<- as.numeric(modelo.actual$fitted.values)
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  h1("Ejemplo para la eleccion de modelos usando LRT"),

    mainPanel(
      tabsetPanel(
        id = 'dataset',
        tabPanel("Exploracion",
                 selectInput('x', 'Variables explicativas X', names(mtcars)[c(1,3, 4,5,6,7)])
                   ,plotlyOutput("plot")),
        
        tabPanel("Seleccion de modelo lineal ",
                 tabItem(tabName="foo",
                         withMathJax(includeMarkdown("Tehory.Rmd"))),
                  hr(),
                 h3('La densidad de las variables independientes es la siguiente'),
                 plotlyOutput("kernel"),
                 h4('El modelo estimado por MV se ve de la siguiente forma'),
                 plotlyOutput("modelo.actual") ,
                 h5("Ejercicio"),
                 h3("Construye modelo lineales con las siguientes variables,"),
                 h3("reporta cuales de ellos aceptarias la hipotesis nula"),
                 h1("Analisis","LRT!"),
                 tableOutput("LRT"),
                 checkboxGroupInput("variables", label = h3("Variables en el modelo"), 
                                    choices = list("disp",  "hp", "drat", "wt", "qsec" ),
                                    selected = "disp")
                ),
        tabPanel("Doc",
                 tabItem(tabName="foo",
                         withMathJax(includeMarkdown("Datos.Rmd"))),
                 hr()),
        tabPanel(
          "Rcpp",
          sidebarLayout(
            sidebarPanel(
              
            ),
            mainPanel(
              
              #uiOutput("testing")
            )
          ) )
        
        
        
      )
    )
  )



# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  
  
   output$plot <-  renderPlotly({ 
     plot_ly(data = mtcars,  x = ~get(input$x), y = ~mpg,
                marker = list(size = 10,
                              color = 'rgba(30,0,125,.6)',#'rgba(30,0,125,.6)',#'rgba(30,0,125,.6)',
                              line = list(color = 'rgb(58,200,225)',
                                          width = 1)),
             text = ~paste("mpg: ", mpg, paste0('<br> ', (input$x),":") ,get(input$x) )) %>%
     layout(title = 'Dispersion de las variables contra mpg')
       
   })
  
     output$modelo.actual <-  renderPlotly({
       graph_reso <- 1
       #Setup Axis
       axis_x <- seq(min(carros$disp), max(carros$disp), by = graph_reso)
       axis_y <- seq(min(carros$drat), max(carros$drat), by = graph_reso)
       #Sample points
       petal_lm_surface <- expand.grid(disp = axis_x, drat = axis_y,KEEP.OUT.ATTRS = FALSE)
       petal_lm_surface$z_hat <- predict(modelo.actual, newdata = petal_lm_surface)
       petal_lm_surface <- acast(petal_lm_surface,  drat~disp,value.var = "z_hat") #y ~ x
       
       output$modelo.actual <-  renderPlotly({ 
         d <- plot_ly(carros, 
                    x = ~disp, 
                    y = ~drat, 
                    z = ~mpg,
                    type = "scatter3d", 
                    mode = "markers",
                    marker = list(size = 6,
                                  color = 'rgba(195,47,47,.9)',#'rgba(30,0,125,.6)',
                                  line = list(color = 'rgb(58,200,225)',
                                              width = .2))) %>% add_trace(
                                z = petal_lm_surface,
                                x = axis_x,
                                y = axis_y,
                                type = "surface") 
       })
       output$kernel <-  renderPlotly({ 
         kd <- kde2d(carros[,"disp"], carros[,"drat"], n = 32)
          plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% 
           add_surface(x = kd$x, y = kd$y, z = kd$z, opacity = 0.98)
         
       })
       
       selectedData.Exploracion <- reactive({
         carro2<- mtcars[, c(input$variables,"mpg" )]
         modelo.reducido <- lm(mpg~ ., carro2)
         m <- lrtest(modelo.reducido, modelo.actual)
         return(m) 
       })
       output$LRT <- renderTable({
         (selectedData.Exploracion())
       })
       # output$testing <- renderUI({
       #   
       #   list(
       #     # Added parameter to `tes1` to pass in krit.
       #     #renderPrint(tes1(n = 3, krit = krit))
       #     
       #     renderPrint(microbenchmark(seq = (em3(y = mY,
       #                                           X = mX,
       #                                           maxit = 100
       #     )
       #     )
       #     ,
       #     para = (em4(y = mY,
       #                 X = mX,
       #                 maxit = 100,
       #                 nthr = 4
       #     )
       #     ),
       #     times = 3
       #     )
       #     )
       #   )
       # })
     })
}



shinyApp(ui, server)