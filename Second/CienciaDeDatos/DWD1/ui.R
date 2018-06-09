####################################
#####    J. Antonio Garcia #########
#####   jose.ramirez@cimat.mx ######
####################################
library(shinydashboard)
library(shiny)
library(plotly)
library(knitr)
library(rmarkdown) 
#########################################
# Construccion de la UI                 # 
#########################################
sidebar <- dashboardSidebar(
  #comenzamos con el menu
  sidebarMenu(
    menuItem("CIMAT", tabName = "CIMAT", icon = icon("dashboard")),
    menuItem("DWD", icon = icon("th"), tabName = "DWD",
             badgeLabel = "nuevo", badgeColor = "green"),
    menuItem("Optimizacion", icon = icon("th"), tabName = "Optimizacion"),
    menuItem("ResultadosDatosSimulados",  icon = icon("th"), tabName = "ResultadosDatosSimulados"),
    menuItem("ResultadosDatosReales",  icon = icon("th"), tabName = "ResultadosDatosReales")
  )
)
#cramos varias tabs
body <- dashboardBody(
  tabItems(
    tabItem(tabName = "CIMAT",
            h2('Distance Weighted Discrimination (DWD)'),
            h3('y la maldición de la dimensionalidad')
    ),
    
    tabItem(tabName = "DWD",
            h2("Geometría de DWD, data piling y el fracaso de SVM"), 
      fluidRow(
            box(  #title = "Distribución sobre la dirección óptima de Bayes",
               #   background = "light-blue",
              #background = "green",
              #solidHeader = TRUE,
              plotlyOutput("puntos", height = 230)
            ),
            box(
             # title = "Distribución sobre la dirección óptima de Bayes",
              #background = "light-blue", solidHeader = TRUE,
              plotlyOutput("Bayes", height = 230)
            )
          ),
            fluidRow(
            box(
              #title = "Proyección sobre la dirección MDP",
               #solidHeader = TRUE, background = 'light-blue',
              plotlyOutput("DWD", height = 230) 
              ),
            fluidRow(
            box(
              title = "¿Qué dimensión?", 
              "Tamaño de muestra fijo 20", br(),
              sliderInput("d", "d", min = 2, max = 1000, step = 5, value = 2)
            ))
          ),
      fluidRow(
          box(width = 12,h4('Las config. HDLSS tienden asintóticamente (d al inf. y n fijo) a tener una estructura geométrica fundamentalmente rígida'), 
            h4('La principal fortaleza de DWD es que su desempeño es cercano al de SVM, cuando SVM es mejor'),
               h5('Cuando d >>n los datos consisten en un subespacio n dimensional y la idea de trabajar en este espacio es impráctica'),
                h6('Los nuevos datos se espera que aparezcan fuera de este subespacio'),
            h6('En el contexto de microarreglos el interés recae en solo algunos subconjuntos de genes específicos y esta atención es más difícil de mantener solo con algunas combinaciones lineales (es decir cualquier base del subespacio generado por los datos) de lo genes considerados')
          ))
    ),
    #la tab de la derivacion
     tabItem(tabName = "Optimizacion",  h2("Problema de optimización de DWD"),
             fluidRow( h1('                        '),
               box(width = 12,  column(6,img(src='margen.png', align = "center", height = 400)),
                       column(4, img(src='cono.png', align = "center", height = 400))
               )), hr(),
             fluidRow( 
               box( width = 12,       column(6,  withMathJax(includeMarkdown(("SVM.Rmd")))) ,
                   
                       column(6, withMathJax(includeMarkdown(("planteamientoDWD.Rmd"))))
                 )
             )
            ),
    #la tab de resultados experimentales
    tabItem(tabName = "ResultadosDatosSimulados",  h1("ResultadosDatosSimulados"),
            fluidRow( h1('                        '),
                      box( column(6, img(src='simulacion1.png', align = "center", height = 500))),
                      box(column(6,  img(src='simulacion2.png', align = "center", height = 500)))), 
            fluidRow(box(width = 12,column(4,
                             img(src='simulacion3.png', align = "center", height = 500)
                      )
                      ))),
    tabItem(tabName = "ResultadosDatosReales",  h1("ResultadosDatosReales"),
            fluidRow( h1('                        '),
                      box( column(5, img(src='cancergrande.jpeg', align = "center", height = 700))),
                      box(column(5,  img(src='cancer_coom.jpeg', align = "center", height = 700)))), 
            fluidRow(box( column(5, img(src='genes_todos.jpeg', align = "center", height = 700))),
                     box(column(5,  img(src='genes_intervalos.jpeg', align = "center", height = 700)))))
            ))



# Put them together into a dashboardPage
dashboardPage(skin = "purple", 
  dashboardHeader(title = "CIMAT Monterrey"),
  sidebar,
  body
)