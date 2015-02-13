library(shiny)
library(MASS)
#library(car)

shinyServer(function(input, output, session) {
  
  s=''
  
  output$distPlot <- renderPlot({
    
    inFile <- input$file1
    
    if (is.null(inFile)){
      #return(NULL)
      #debug
      m=read.csv("rt.csv")
      x=m$v
    }else{
      m=read.csv(inFile$datapath)
      x=as.matrix(m[-1])  #m$v
    }
    
    xMin=0
    xMax=input$valid[2] #max(1,max(x))
    bins <- seq(min(x), xMax, by=input$tstep)
    
    # draw the histogram with the specified number of bins
    h=hist(x, breaks = bins, plot=T, probability=TRUE, main=sprintf('%s',input$model))
    #h$counts=h$counts/sum(h$counts)
    #plot(h, xlim = c(0,xMax), 
    #     col = 'gray', border = 'darkgreen')
  
    
    #model
    x[x<input$valid[1]]=NA
    x[x>input$valid[2]]=NA
    
    
    abline(v=quantile(x, c(0.5)), col='blue')
    abline(v=quantile(x, c(0.25, 0.75)), col='magenta')
    M=mean(x)
    abline(v=M, col='green')
    
    xx=seq(xMin, xMax, by=0.004)
    
    #text(x=0,y=0.1,labels = input$model)
    
    
    #method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), 
    s=''
    if (is.element(input$model, c('Density'))){
      fn=(density(x, bw=input$bw, na.rm=T))
      xx=fn$x
      yy=fn$y
      s= sprintf('bandwidth = %.3f, kernel="%s"', fn$bw,'gaussian')
    
    
    }else if (input$model=='Weibull-Gaussian'){
      
      #fn= function(x,m,s,v,l) (l/v)*(x/v)^(l-1) * exp(-x/v) * 1/((2*pi)^0.5 * s) * exp(-(x-m)^2/(2*s^2))
      #p=fitdistr(x, fn, start=list(m=0.3, s=0.1, v=0.1, l=3.0), lower=0)
      
      fn= function(x, m,s,b,c)  c * b^(-c) * x^(c-1) * exp(-(x/b)^c)  * 1/((2*pi)^0.5 * s) * exp(-(x-m)^2/(2*s^2))
      p=fitdistr(x, fn, start=list(m=0.3, s=0.1, b=0.25, c=5.2), 
                 control = list(trace = 1,fnscale=-1),
                 lower=c(0.1, 0.01, 0.12, 1.1),upper = c(1.0,1.0,1,10))
      
      kk=round(coef(p), 4)
      yy=fn(xx, kk[1],kk[2],kk[3],kk[4])
      
      s=expression(paste(mu%==% sigma, kk[1],kk[2],kk[3],kk[4]))
      
    }else if (input$model=='ex-Gaussian'){
        #fn=function(t,m,s,v){
          #o=integrate(function(y) exp(-y/v)/(v) * 1/((2*pi)^0.5 * s) * exp(-(t-y-m)^2/(2*s^2))), 0, Inf)
          #return(o$value)
        #}
        #-log(nu)+((mu-x)/nu)+(sigma^2/(2*nu^2))+log(pnorm(((x-mu)/sigma)-sigma/nu))
        fn=function(t,mu,sigma,nu){
          z=t-mu-((sigma^2)/nu)
          return( exp(-log(nu)-(z+(sigma^2/(2*nu)))/nu+log(pnorm(z/sigma))))
        }
        p=fitdistr(x, fn, start=list(mu=0.3, sigma=0.03, nu=0.1), lower=c(0.1,0.01,0.03), upper=c(xMax, 1, 2))
        kk=round(coef(p), 3)
        yy=fn(xx, kk[1], kk[2], kk[3])
        #yy=yy/sum(yy)
        lines(xx, yy, col='black')
        
        s=(sprintf('mu = %.3f, sigma = %.03f, nu = %.3f ',
                   kk[1], kk[2], kk[3])
        )
        
        
        #when nu<0.05
        #dnorm(x, mean=mu, sd=sigma, log=TRUE)
        
        #exgaussalt
        #f(x; y_0, A, x_c, w, t_0 )=y_0+\frac{A}{t_0} \exp \left( \frac {1}{2} \left( \frac {w}{t_0} \right)^2 - \frac {x-x_c}{t_0} \right) \left( \frac{1}{2} + \frac{1}{2} \operatorname{erf} \left( \frac {z}{\sqrt{2}} \right) \right) ,
        
        
    } else{
      p=fitdistr(x, tolower(input$model))
      kk=round(coef(p), 3)
      yy=switch(input$model,
          #fn = function(x) 1/(gamma(a)) * b^a * x^(a-1) * exp(-xb)
          Gamma=dgamma(xx, p$estimate[1], p$estimate[2]),
          Normal=dnorm(xx, p$estimate[1], p$estimate[2]),
          Weibull=dweibull(xx, p$estimate[1], p$estimate[2])
          )
      s=sprintf('parameters: %.3f,  %.3f', kk[1], kk[2])
    }

## goodness of fit???
#check http://cran.r-project.org/doc/contrib/Ricci-distributions-en.pdf    
#qqPlot(x, distribution="weibull", shape=5.2, scale=0.327)

    #title(expression(s))
    #h$counts=h$counts/sum(h$counts)
    #plot(h, xlim = c(0,xMax), 
    #     col = 'gray', border = 'darkgreen')
    lines(xx, yy, col='black', lwd=2)
    
    mtext(s, side=3, adj=1)
    
  })
  

  

  output$info <- renderUI({
    
    HTML(switch(input$model,

Density='Density approximation. You can adjust band width with control.',
Normal='<p>Symmetric Gaussian distribution around (1) mean with (2) variation between 
      -3&sigma; and 3&sigma;</p>
<p>
          $$f(x)=\\frac{1}{\\sqrt{2}\\sigma} \\int e^{-t} dt$$
</p>',
          
Gamma='Asymmetric distribution. Parameters: (1) shape, (2) scale
          ',

Weibull='Asymmetric distribution. Parameters: (1) shape, (2) scale

        $$ f(x) = \\frac{c}{b^c} x^{c-1} e^{-\\left(\\frac{x}{b}\\right)^c} $$

          ',

"ex-Gaussian"='ex-Gaussian
    $$f(x)=\\frac{\\lambda}{\\nu}\\left(\\frac{\\sigma^2}{2*\\nu}\\right)^{\\lambda-1}$$
    <p>
    See also:<br/>
          Sternberg S. Reaction Times and the Ex-Gaussian Distribution: When is it Appropriate?<br>
          Sternberg, S. (2014). Sequential processes and the shapes of reaction-time distributions.
    </p>
',

'Weibull-Gaussian'='<p>The Weibull-Gaussian is a convolution of the Gaussian (normal) and Weibull.</p> 
    $$f(x)=\\frac{\\lambda}{\\nu}\\left(\\frac{y}{\\nu}\\right)^{\\lambda-1} + \\left(\\frac{\\sigma^2}{2*\\nu}\\right)^{\\lambda-1}$$
<img "src=Palmer_Horowitz_2011.png"/>
    '
          ), '<script>MathJax.Hub.Queue(["Typeset", MathJax.Hub]);</script>')
  })
  
})

