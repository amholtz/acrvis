library(tidyverse, quietly = T, verbose = F, warn.conflicts = F)
library(sf, quietly = T, verbose = F)
library(tidytree, warn.conflicts = F, quietly = T)
library(shiny)
library(leaflet)
source(here::here("map_migrationRates.R"))



options(shiny.maxRequestSize = 30*1024^2)

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Ancestral Character Map Visualization"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select log file ----
      fileInput("file1", "Choose log File",
                multiple = FALSE
      ),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Select gps file ----
      fileInput("file2", "Choose GPS File",
                multiple = FALSE
      ),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Select mcc file ----
      fileInput("file3", "Choose MCC File",
                multiple = FALSE
      )
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      # h3("Log File Sample"),
      #  tableOutput("contents_log"),
      #  h3("GPS File Sample"),
      #  tableOutput("contents_gps"),
      #  h3("MCC File Sample"),
      #plotOutput("contents_mcc"),
      h5("Plot migration rates with BF >= 3"),
      leafletOutput("migration_rates"),
      
    ))
)






# Define server logic to read selected file ----
server <- function(input, output) {
  
  
  
  #LOG RENDERING
  output$contents_log <- renderTable({
    
    req(input$file1)
    tryCatch(
      {
        df <- read.table(input$file1$datapath,
                         header = T,
                         sep = "\t")
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    return(head(df))
    
  })
  
  #GPS RENDERING
  output$contents_gps <- renderTable({
    
    req(input$file2)
    tryCatch(
      {
        df <- read.delim(input$file2$datapath, header=FALSE)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    return(head(df))
    
  })
  
  #MCC RENDERING
  #output$contents_mcc <- renderPlot({
    
  #  req(input$file3)
   # tryCatch(
   #   {
   #     mcc <- read.beast(input$file3$datapath)
   #   },
   #   error = function(e) {
        # return a safeError if a parsing error occurs
   #     stop(safeError(e))
   #   }
   # )
    
   # return(ggtree(mcc))
    
 # })
  
  
  output$migration_rates <- renderLeaflet({
    
    req(input$file1)
    req(input$file2)
    
    log <- read.table(input$file1$datapath,
                      header = TRUE,
                      sep = "\t")
    gps <- read.delim(input$file2$datapath, header=FALSE)
    
    maps_migrationrates(log = log, gps = gps)
    
  })
}

# Create Shiny app ----
shinyApp(ui, server)