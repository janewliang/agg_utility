library(shiny)
library(shinybrowser)
library(rhandsontable)

library(plotly)

# Define UI for application
shinyUI(fluidPage(
    detect(),

    # Application title
    titlePanel("Aggregate utility"),
    br(), 
    
    fluidRow(
        # Instructions
        column(5, 
               p("Enter the prevalence, penetrance, and precision estimates for 
                 up to 10 genes in the table."), 
               tags$ul(
                   tags$li("Right-click the table to add/remove rows."), 
                   tags$li("Empty rows and rows with missing information will 
                           be ignored."), 
                   tags$li("Invalid table values will raise an error."), 
                   tags$li("The test utility (regardless of results) is assumed 
                           to be 0.")),
               p("Hover over the plots to see values for the net utility 
                 (estimate, 2.5%, and 97.5%) and the probability of the net 
                 utility being positive at different different FN/FP utility 
                 cost ratios.")), 
        
        # Table
        column(7, rHandsontableOutput("table_of_values"), 
               tags$head( # Remove horizontal scroll
                   tags$style(HTML(".handsontable {
                                   overflow: visible;
                                   }"))
               ))), 
    
    fluidRow(
        # Link to Github repo
        column(12, 
               "Source code: ", 
               tags$a(href="https://github.com/janewliang/agg_utility", 
                      "github.com/janewliang/agg_utility"))
    ), 
    br(), 
    
    fluidRow(
        # Plots
        column(12, align="center", plotlyOutput("plots"))
        )
))
