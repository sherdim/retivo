library(shiny)
library(MASS)

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
    
    tStep=0.004
    xx=seq(xMin, xMax, by=tStep)
    
    #text(x=0,y=0.1,labels = input$model)
    
    
    #method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), 
    s=''
    if (is.element(input$model, c('Density'))){
      fn=(density(x, bw=input$bw, na.rm=T))
      xx=fn$x
      yy=fn$y
      xM=xx[which.max(yy)]
      s= sprintf('bandwidth = %.3f, kernel="%s", mode = %.3f', fn$bw,'gaussian',xM)
    
    
    }else if (input$model=='Weibull-Gaussian'){
      #http://stackoverflow.com/questions/23569133/adding-two-random-variables-via-convolution-in-r 
      
      N = function(x,m,s) dnorm(x,m,s)        # normal
      W = function(x,v,l) dweibull(x,v,l)   # weibull( shape, scale
      Z=function(x,m=0.3,s=0.1,v=1.2,l=1) W(x,v,l) * N(x,m,s)
      fn=Vectorize(Z)
      
      
      #op <- options(digits = 3)
      #fn= function(x,m,s,v,l) (l/v)*(x/v)^(l-1) * exp(-x/v) * 1/((2*pi)^0.5 * s) * exp(-(x-m)^2/(2*s^2))
      p=fitdistr(x, fn, start=list(m=0.3, s=0.05, v=5, l=.3), 
                 #control = list(trace = 1,fnscale=1),
                 lower=c(0.1, 0.01, 1.01, 0.1),upper = c(1.0,.3,20,0.9))
      
      kk=coef(p)
      yy=fn(xx, kk[1],kk[2],kk[3],kk[4])
      yy=yy*length(yy)/sum(yy)
      s=sprintf('mu=%.3f, sigma=%.3f, k=%.3f, lambda=%.3f', kk[1],kk[2],kk[3],kk[4])
  

    }else if (input$model=='Zaitsev-Skorik'){
        fn=function(x,a=0.02,b=0.3,k=1) k/(a) * exp(-exp((b-x)/a)+(b-x)/a)
        #p=fitdistr(x, fn, start=list(a=0.02,b=0.3), lower=c(0.001,0.1), upper=c(1,xMax))
        #kk=round(coef(p), 3)
        
        y=ecdf(x)
        #y0=log(log(1.0/y))
        fcdf = function(p, x) exp(-exp((p[2]-x)/p[1]))
        fmin = function(p, x) sum((fcdf(p,x) - y(x))^2)
        fit=optim(list(a=0.02,b=0.3), fmin, x=x)
        kk=fit$par
        yy=fn(xx, kk[1], kk[2]) #/ input$tstep
        #yy=yy/sum(yy)
        lines(xx, yy, col='black')
        
        s=sprintf('a=%f, b=%f', kk[1],kk[2])
        
        
    }else if (input$model=='Wald'){
        fn = function(x, mu=0.3, l=1) {
          return( sqrt(l/2*pi*x^3) * exp(-l*(x-mu)^2/(2*x*mu^2)))
        }
        p=fitdistr(x, fn, start=list(mu=0.3, l=1), lower=c(0.1,0.1), upper=c(xMax, 1000))
        kk=round(coef(p), 3)
        
        yy=fn(xx, kk[1], kk[2])
        yy=yy*length(yy)/sum(yy)
        lines(xx, yy, col='black')
        
        s=sprintf('mu=%f, lambda=%f', kk[1],kk[2])
        
        
    }else if (input$model=='ex-Wald'){
        W = function(x, mu=0.3, l=1) {
          return( sqrt(l/2*pi*x^3) * exp(-l*(x-mu)^2/(2*x*mu^2)))
        }
        X =function (x, lambda=1) dexp(x,lambda)
        
        Z=function(x,mu=0.3,l=1,lambda=1) X(x,lambda) * W(x,mu,l)
        fn=Vectorize(Z)
        p=fitdistr(x, fn, start=list(mu=0.3, l=0.1, lambda=1), lower=c(0.1,0.001,0.0001), upper=c(xMax, 20, 1000))
        
#         fn=function(x,mu=0.3,sigma=0.1,a=1,gamma=0.001){
#           k=sqrt(mu^2 - 2* gamma * sigma^2)
#           #return( gamma * exp(- gamma*t+(a*(mu-k))/sigma^2)) # * pnorm((x-mu)/sigma))
#           return( exp( log(gamma) + (- gamma*t+(a*(mu-k))/sigma^2) + log(pnorm((x-mu)/sigma))))
#         }
        #p=fitdistr(x, fn, start=list(mu=0.3, sigma=0.1, a=1, gamma=0.001), lower=c(0.1,0.001,0.3,0.0001), upper=c(xMax, 1, 1000, 1))

        kk=coef(p)
        yy=fn(xx, kk[1], kk[2], kk[3])
        yy=yy*length(yy)/sum(yy)
        lines(xx, yy, col='black')
              
        s=sprintf('mu = %.3f, lambda = %.03f, gamma = %.3f ', kk[1], kk[2], kk[3])


    }else if (input$model=='ex-Gaussian'){
        #fn=function(t,m,s,v){
          #o=integrate(function(y) exp(-y/v)/(v) * 1/((2*pi)^0.5 * s) * exp(-(t-y-m)^2/(2*s^2))), 0, Inf)
          #return(o$value)
        #}
        #-log(nu)+((mu-x)/nu)+(sigma^2/(2*nu^2))+log(pnorm(((x-mu)/sigma)-sigma/nu))
        fn=function(t,mu=0.3,sigma=0.03,nu=0.1){
          z=t-mu-((sigma^2)/nu)
          return( exp(-log(nu)-(z+(sigma^2/(2*nu)))/nu+log(pnorm(z/sigma))))
        }
        p=fitdistr(x, fn, start=list(mu=0.3, sigma=0.03, nu=0.1), lower=c(0.1,0.01,0.03), upper=c(xMax, 1, 2))
        kk=round(coef(p), 3)
        yy=fn(xx, kk[1], kk[2], kk[3])
        #yy=yy/sum(yy)
        lines(xx, yy, col='black')
        
        s=(sprintf('mu = %.3f, sigma = %.03f, nu = 1/lambda = %.3f ',
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

$$f(t)=\\frac{(x-\\xi)^{\\beta-1} e^{-\\frac{x-\\xi}{a}}}{a^{\\beta} \\Gamma(\\beta)}$$
          ',

Weibull='Asymmetric distribution. Parameters: (1) shape, (2) scale

        $$ f(t) = \\frac{c}{b^c} t^{c-1} e^{-\\left(\\frac{t}{b}\\right)^c} $$

          ',

"ex-Gaussian"='ex-Gaussian - convolution of a Gaussian and an exponential distribution,
named by Burbeck & Luce (1982)

    $$f(x)=\\frac{1}{\\lambda} e ^{\\frac{\\mu}{\\lambda}+\\frac{\\sigma^2}{2*\\lambda^2}-\\frac{x}{\\lambda}} \\Phi \\left(\\frac{x-\\mu-\\sigma^2/\\lambda}{\\sigma}\\right)$$

    <p>
    See also:<br/>
<li>
Burbeck, S. L., & Luce, R. D. (1982). Evidence from auditory simple
reaction times for both change and level detectors.
Perception & Psychophysics, 32 (2), 117-133.
<li>Sternberg S. (2014) Reaction Times and the Ex-Gaussian Distribution: When is it Appropriate?<br>
<li>Sternberg, S. (2014). Sequential processes and the shapes of reaction-time distributions.
<li>Van Zandt, T. (2000). How to fit a response time distribution.
Psychonomic Bulletin &amp; Review, 7 (3), 424-465.
<li>Hwang Gu SL, Gau SS, Tzang SW, Hsu WY. The ex-Gaussian distribution of reaction times in adolescents with attention-deficit/hyperactivity disorder. Research in Developmental Disabilities. 2013 Nov;34(11):3709-3719. 
    </p>
',

'Weibull-Gaussian'='<p>The Weibull-Gaussian is a convolution of the Gaussian (normal) and Weibull distributions.</p> 

$$f(t)=\\int \\frac{\\lambda}{\\nu}\\left(\\frac{y}{\\nu}\\right)^{\\lambda-1} e^{\\frac{y}{\\nu}} \\frac{1}{\\sqrt{2*\\pi}\\sigma}  e^{-\\frac{(t-y-\\mu)^2}{2\\sigma^2}} dy$$

<p>
<img src="Matsushita_2007.png"/></p>

<li>Matsushita, S. (2007). RT analysis with the Weibull-Gaussian convolution model. Proceedings of Fechner Day, 23(1).
<li>Marmolejo-Ramos, F., & González-Burgos, J. (2013). A power comparison of various tests of univariate normality on Ex-Gaussian distributions. Methodology: European journal of research methods for the behavioral and Social sciences, 9(4), 137.

',

'Wald'=' inverse Gaussian distribution (also known as the Wald distribution) 

$$f(x)= \\sqrt{\\frac{\\lambda}{2 \\pi x^3}} e^{\\frac{-\\lambda(x-\\mu)^2}{2x\\mu^2}} $$
',

'ex-Wald'='ex-Wald is the convolution of an exponential and a Wald
distribution that attempts to represent both the decision and response components of a response time (Schwarz, 2001) as a
diffusion process, a continuous accumulation of information towards some decision threshold

$$f(t) =  \\gamma exp\\left(-\\gamma t + \\frac{a(\\mu-k)}{\\sigma^2}\\right) * F(t)$$

$$ k = \\sqrt{\\mu^2 - 2\\gamma\\sigma^2} \\geq 0 $$

<p>
<img src="Palmer_Horowitz_2011.png"/></p>

<li>Palmer, E. M., Horowitz, T. S., Torralba, A., & Wolfe, J. M. (2011). What are the shapes of response time distributions in visual search? Journal of Experimental Psychology: Human Perception and Performance, 37.
<li>Schwarz, W. (2001). The ex-Wald distribution as a descriptive model of
response times. Behavioral Research Methods, Instruments & Computers, 33 (4), 457-469.

    ',

'Zaitsev-Skorik'='<p>
$$f(t)=\\frac{\\lambda}{a} exp\\left[-exp(\\frac{b-x}{a}) + \\frac{b-x}{a}\\right]$$
</p>
<p>
a - razbros, b - moda
</p>
<li>Zaitsev, Skorik, 2002 Mathematical Description of Sensorimotor
Reaction Time Distribution, Human Physiology, No.4

'
          ), '<script>MathJax.Hub.Queue(["Typeset", MathJax.Hub]);</script>')
  })
  
})

