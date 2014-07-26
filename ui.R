library(shiny)

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("RT distribution analysis"),
  
  # Sidebar with a slider input for number of bins
  sidebarPanel(
    fileInput('file1', 'Choose file to upload',
              accept = c(
                'text/csv',
                'text/comma-separated-values',
                'text/tab-separated-values',
                'text/plain',
                '.csv',
                '.tsv'
              )
    ),
    p('Pick .csv or .tsv file with two columns: time and reaction time, s.' ),
    tags$hr(),
    sliderInput("valid",
                "Valid diapason, s:",
                min = 0.0,
                max = 3.0,
                #round=-3,
                value = c(0.1, 1)),
    sliderInput("tstep",
                "Size of bins, s:",
                min = 0.01,
                max = 0.1,
                step = 0.005,
                value = 0.025),
    
    tags$hr(),
    p('Parameters of distribution model'),
    selectInput('model','Model Type:',
    c(density='density','weibull',ecdf='Empirical CDF','Weibull-Gaussian')),
    sliderInput("skewness",
                "Skewness:",
                min = 0,
                max = 1,
                value = 0.2)
    
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    plotOutput("distPlot")
  )
))
