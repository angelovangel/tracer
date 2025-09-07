library(shiny)
library(bslib)
library(dplyr)
library(reactable)
library(bsicons)
library(purrr)
library(sparkline)
#library(sangerseqR)

source('global.R')

sidebar <- sidebar(
    fileInput("ab1", "Choose ab1 files", accept = ".ab1", multiple = T),
)

ui <- page_navbar(
  sidebar = sidebar,
  title = "Tracer - Sanger reads viewer",
  #nav_spacer(),
  #nav_item(selectInput('test', '', c('a', 'b'), multiple = T)), 
  nav_panel('',
    reactableOutput('table')
  )
)

server <- function(input, output, session) {
  
  
  df <- reactive({
    ab1list <- input$ab1$datapath
    req(ab1list)
    purrr::map_dfr(ab1list, get_ab1)
  }) 
    
  
  output$table <- renderReactable({
    data <- df() %>% 
      select('sample', 'rundate', 'rawSeqLen', 'trimSeqLen', 'crl20', 'seq', 'seqtrimmed')
    reactable(data)
  })
}

shinyApp(ui, server)