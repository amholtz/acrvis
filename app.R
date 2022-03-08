
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.14")


if (!"treeio" %in% installed.packages()) {
  if (!"devtools" %in% installed.packages()) install.packages("devtools")
  devtools::install_github("GuangchuangYu/treeio")
} else if (packageVersion("treeio") < "1.5.1.2") {
  if (!"devtools" %in% installed.packages()) install.packages("devtools")
  devtools::install_github("GuangchuangYu/treeio")
}

if (!"purrr" %in% installed.packages()) install.packages("purrr")


purrr::walk(c("shiny", "shinyjs", "tidyverse", 
              "tidytree", "shinyalert"), ~{
                if (!.x %in% installed.packages()) install.packages(.x)
              })


suppressWarnings({
  suppressPackageStartupMessages({
    require(shiny)
    require(sf)
    require(shinyjs)
    require(tidyverse)
    #require(ggtree)
    require(tidytree)
    require(treeio)
    require(shinyalert)
    require(leaflet)
  })
})

require(BiocManager)
options(repos = BiocManager::repositories())

source(here::here("map_migrationRates.R"))
#source(here::here("read_beast.R"))


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
      h3("MCC File Sample"),
      plotOutput("contents_mcc"),
      h5("Plot migration rates with BF >= 3"),
      leafletOutput("migration_rates")
      
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
  output$contents_mcc <- renderPlot({
    
    req(input$file3)
    tryCatch(
      {
        mcc <- treeio::read.beast(input$file3$datapath)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    
    return(ggtree::ggtree(mcc))
    
  })
  
  
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