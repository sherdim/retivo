library(shiny)

zz=c('Normal', 'Weibull','Gamma',
     'Density',
     'Wald',
     'ex-Gaussian',#='exgaussian',
     'ex-Wald',
     'Weibull-Gaussian',#='weibullgaussian'
   'Zaitsev-Skorik'
)

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
    p('Pick .csv or .tsv file with two columns: time and reaction time, s.',
      a('example',href='rt1.csv')),
    tags$hr(),
    sliderInput("valid",
                "Valid diapason, s:",
                min = 0.0,
                max = 3.0,
                step = 0.005,
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
    selectInput('model','Model Type:', zz),
    conditionalPanel(condition = "input.model == 'Density'",  
      sliderInput("bw",
                  "Bandwidth (density):",
                  min = 0.005,
                  max = 0.1,
                  step = 0.005,
                  value = 0.02)
    #),
    #conditionalPanel(condition = "input.model == 'Normal'",  
      #sliderInput("skewness",
       #         "Skewness:",
        #        min = -1,
         #       max = 1,
          #      value = 0)
    )
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    plotOutput("distPlot"),
    withMathJax(p('Green line - mean, Blue line - median, Magenta lines -  25 and 75 quantiles')),
    #uiOutput('formula'),
    
    #h3('Information about selected distribution model'),
    uiOutput('info')
    #verbatimTextOutput('info')
  )
))
