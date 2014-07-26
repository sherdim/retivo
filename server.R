library(shiny)

shinyServer(function(input, output) {
   
  output$distPlot <- renderPlot({
    
    inFile <- input$file1
    
    if (is.null(inFile)){
      #return(NULL)
      #debug
      m=read.csv("~\\retivo\\rt.csv")
      x=m$v
    }else{
      m=read.csv(inFile$datapath)
      x=as.matrix(m[2])  #m$v
    }
    
    xMin=0
    xMax=input$valid[2] #max(1,max(x))
    bins <- seq(min(x), xMax, by=input$tstep)
    
    # draw the histogram with the specified number of bins
    h=hist(x, breaks = bins, plot=F)
    h$counts=h$counts/sum(h$counts)
    plot(h, xlim = c(0,xMax), 
         col = 'gray', border = 'darkgreen')
    
    #model
    x[x<input$valid[1]]=NA
    x[x>input$valid[2]]=NA
    
    
    
    text(x=0,y=0.2,labels = input$model)
    z=switch(input$model,
          density=density(x),
          ecdf=ecdf(x)
           )
    #Y=Weibull-Gaussian
    
    text(x=0,y=0.15,labels = input$skewness)
    
    #X=seq(xMin, xMax, 0.005)
    #Y=predict(z, X)
    #plot(z, col='black')
    
    s='$y=exp(x)$'
    text(x=0.5,y=0.2,labels = s)
    
    text(x=0.5,y=0.1,labels = z)
    
  })
  
  
  output$m <- renderTable({
    
  })
  
})

