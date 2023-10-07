# Two Way HANOM -----------------------------------------------------------

#install.packages("pacman") #Only first time needed
#devtools::install_github('Mikata-Project/ggthemr') #Only first time needed
pacman::p_load(shiny, dplyr, data.table, shinythemes, shinycssloaders, waiter, ggplot2, ggtext, ggthemr, waiter)

# Functions ---------------------------------------------------------------
# Calculate Statistics
MeanP1<- function(data){
    
    #Current Restriction: 1st column as value, 2nd column as treatment/ factor
    data<- data.frame(data)
    tr<- as.factor(data[,2])
    val<- as.numeric(data[,1])
    
    k <- levels(tr)
    nk <- nlevels(tr)    #levels of treatment
    ni <- tr %>% table #number of samples in each treatment
    
    # Random sample ni-1 ro calculate sample mean and sample variance
    xbari <- sapply(1:nk, function(a)  mean(val[tr == k[a]][-ni[a]]))
    vari <- sapply(1:nk, function(a) var(val[tr == k[a]][-ni[a]]))
    
    # Calculate weighted sample mean: Xtilde
    maxvar <- max(vari)
    Ui <- sapply(1:nk, function(a) (1/ni[a]) + (1/ni[a]) * sqrt(max(0, (1/(ni[a]-1))*(maxvar/vari[a] - 1)) ) )#Eq4
    Vi <- sapply(1:nk, function(a) (1/ni[a]) - (1/ni[a]) * sqrt(max(0, (ni[a]-1)*(maxvar/vari[a] - 1)) ) )  #Eq5
    xtildei <- sapply(1:nk, function(a) Ui[a] * sum(val[tr == k[a]][-ni[a]]) + Vi[a] * val[tr == k[a]][ni[a]])
    
    out<- rbind(ni, xbari, vari, Ui, Vi, xtildei)
    return(out)
}
MeanP2<- function(data){
    
    #Current Restriction: 1st column as value, 2nd column as treatment/ factor
    data<- data.frame(data)
    tr<- as.factor(data[,2])
    val<- as.numeric(data[,1])
    
    k <- levels(tr)
    nk <- nlevels(tr)    #levels of treatment
    #k<- unique(tr)
    #nk<- length(k)
    ni <- tr %>% table  #number of samples in each treatment
    n0<- min(ni-1)
    
    # Random sample ni-1 ro calculate sample mean and sample variance
    xbari <- sapply(1:nk, function(a) mean(val[tr==k[a]][c(1:n0)], na.rm=TRUE) )
    vari <- sapply(1:nk, function(a)  var(val[tr == k[a]][c(1:n0)], na.rm=TRUE) )
    
    # Calculate weighted sample mean: Xtilde
    maxvar <- max(vari/ni)
    Ui <- sapply(1:nk, function(a) (1/ni[a]) + (1/ni[a]) * sqrt( max(0, ((ni[a]-n0)/n0) *(ni[a]*maxvar/vari[a] - 1)) )) #Eq15
    Vi <- sapply(1:nk, function(a) (1/ni[a]) - (1/ni[a]) * sqrt( max(0, (n0/(ni[a]-n0) * (ni[a]*maxvar/vari[a] - 1))) )) #Eq16
    xtildei <- sapply(1:nk, function(a) Ui[a] * sum(val[tr == k[a]][c(1:n0)]) + Vi[a] * sum(val[tr == k[a]][c((n0+1):ni[a])]) )
    
    out<- rbind(ni, n0, xbari, vari, Ui, Vi, xtildei)
    return(out)
}

# Distribution functions
hdist<- function(ni, iter){
  k<- length(ni)
  ti<- sapply(1:k, function(p) rt(iter, ni[p]-2))   #t(ni-2) distribution
  ttildei<- sapply(1:iter, function(i) {
    sapply(1:k, function(p) ((k-1)/k)*ti[i,p]-(sqrt(ni[p])/k)*sum(ti[i,-p]/sqrt(ni[-p])))  #Eq9 (transform student-t dis into ttilde)
  })
  cv<- data.frame(Min= apply(ttildei, 2, min), Max= apply(ttildei, 2, max))
  return(cv)
}
ddist<- function(n0, iter){
  k<- length(n0)
  ti <- sapply(1:iter, function(i) rt(n = k, df = n0[1]-1))
  out<- apply(ti, 2, function(t) max(t - mean(t)))
  return(out)
}

# Compute P-Value 
# W must be a WeightedMean() object
# cvdist must be a SimCV() object
pvalP1<- function(W, cvdist){
    xbari<- W['xbari',]
    vari<- W['vari',]
    xtildei<- W['xtildei',]
    ni<- W['ni',]
    
    tval<- (xtildei- mean(xtildei))/(sqrt(max(vari)/ni))
    
    pseq<- seq(0, 1, by= 0.001)
    pquantile<- sapply(pseq, function(w) quantile(cvdist, 1-w/2))
    pequation<- abs(pquantile- max(abs(tval)))
    pval<- pseq[pequation== min(pequation)]
    
    return(pval)

}
pvalP2<- function(W, cvdist){
    xbari<- W['xbari',]
    vari<- W['vari',]
    xtildei<- W['xtildei',]
    ni<- W['ni',]
    
    tval<- (xtildei- mean(xtildei))/(sqrt(max(vari/ni)))
    
    pseq<- seq(0, 1, by= 0.001)
    pquantile<- sapply(pseq, function(w) quantile(cvdist, 1-w/2))
    pequation<- abs(pquantile- max(abs(tval)))
    pval<- pseq[pequation== min(pequation)]
    
    return(pval)
}

MetricsP1<- function(W, criticalvalue){
    xtildei <- W["xtildei",]
    k<- ncol(W)
    ni <- W["ni",]
    z1<- max(W["vari",]) 
    
    center <- mean(xtildei)
    LDL <- mean(xtildei) - criticalvalue * sqrt(z1) / sqrt(ni)
    UDL <- mean(xtildei) + criticalvalue * sqrt(z1)/ sqrt(ni)
    
    out<- rbind(W, center=center, LDL=LDL, UDL=UDL)
    #rownames(out) <- c("\\( n_{i} )\\", "\\( \\bar{X_i} \\)", "\\( S_{i} \\)", "\\( U_{i} \\)", "\\( V_{i} \\)", "\\( \\tilde{X_i} \\)", "Center", "LDL", "UDL")
    
    return(out)
}
MetricsP2<- function(W, criticalvalue){
    xtildei <- W["xtildei",]
    k<- ncol(W)
    ni <- W["ni",]
    z1<- max(W["vari",]/ni)

    center <- mean(xtildei)
    LDL <- mean(xtildei) - criticalvalue * sqrt(z1)
    UDL <- mean(xtildei) + criticalvalue * sqrt(z1)
    
    out<- rbind(W, center=center, LDL=LDL, UDL=UDL)
    #rownames(out) <- c("\\( n_{i} )\\", "\\( \\bar{X_i} \\)", "\\( S_{i} \\)", "\\( U_{i} \\)", "\\( V_{i} \\)", "\\( \\tilde{X_i} \\)", "Center", "LDL", "UDL")
    
    return(out)    
}

ggChart2 <- function(charttable, text1, text2, text3){
    k<- colnames(charttable)
    M<- data.frame(t(charttable))
    M$xlabel <- k
    M$xnum <- as.factor(k) %>% as.numeric
    q<- nrow(M)
    
    g<- ggplot(M, aes(xnum, xtildei, group=1)) +
        geom_point() +
        geom_segment(aes(xend= xnum, yend = pmin(xtildei, center))) +
        geom_segment(aes(xend= xnum, yend = pmax(xtildei, center))) +
        
        #coord_cartesian(xlim=c(0.8,q+0.2)) +
        geom_hline(aes(yintercept= center)) + 
        
        # Draw boundaries
        geom_step(aes(x=xnum-0.5, y= UDL), linetype="dashed", linewidth=0.5) +
        geom_step(aes(x=xnum-0.5, y= LDL), linetype="dashed", linewidth=0.5) +
        geom_segment(aes(x=0, xend =0.5, y= UDL[1], yend= UDL[1]), linetype= "dashed", linewidth=0.5) +
        geom_segment(aes(x=Inf, xend = q-0.5, y= UDL[q], yend= UDL[q]), linetype= "dashed", linewidth=0.5) +
        geom_segment(aes(x=0, xend = 0.5, y= LDL[1], yend= LDL[1]), linetype= "dashed", linewidth=0.5)  +
        geom_segment(aes(x=Inf , xend =q-0.5, y= LDL[q], yend= LDL[q]), linetype= "dashed", linewidth=0.5)
    
    
    # Set x-limits
    g<- g + coord_cartesian(xlim=c(0.8,q+0.2)) + 
        scale_x_continuous(breaks= c(1:q), labels = M$xlabel) +
        labs(title=text1, subtitle= text2, x= text3, y= "Xtilde")+
        theme(plot.margin= margin(20,50,25,25), 
              plot.caption = element_text(hjust=0),
              axis.text.x = element_text(size=14)) 
    ggthemr('pale')
    return(g)
}

interaction_test<-function(data){
    colnames(data)[1]<- "value"
    out_aov<- aov(value~.^2, data)
    out2<- data.frame(summary(out_aov)[[1]])
    row_name<- rownames(out2)
    p_val<- out2[grepl(":", row_name),5]
    inter<- sum(which(p_val<= 0.05)) >0  #Interaction effect is significant if result is TRUE
    return(list(result= summary(out_aov)[[1]], sig= inter))
}





# UI ----------------------------------------------------------------------
ui <- fluidPage(
    #autoWaiter(html= spin_pulsar(), color="black"),
    useWaiter(),
    waiterPreloader(html=spin_square_circle(), color="black"),
    #waiterOnBusy(html= spin_dots()),
    #useHostess(),
    withMathJax(),
    theme = shinytheme("united"),
    h1(id = "Two-Way HANOM",   # Application title
       HTML('<span style="color:Tomato; font-family: Optima; font-size: 24;">Two-Way HANOM</span>')
    ),  

    sidebarPanel(
        fileInput("upload","Upload data", multiple = FALSE),
        helpText("1. Keep treatment label names short. Output charts and tables may become cluttered with lengthy label names or numerous groups."),
        helpText("2. Continuous variable must be placed at the first column:"),
        tableOutput("sampledata"),
        br(),
        br(),
        actionButton("run","Run Analysis", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
        width=3,
        br(), br(), br(), br()
        ),
 
    mainPanel(
        numericInput("alpha","Type I error rate (\\(\\alpha\\)):", value= 0.05, min=0, max=1, step= 0.01),
        tabsetPanel(
            id="tabs",
            tabPanel(title="P1 Procedure", 
                     fluidRow(
                         column(8, withSpinner(plotOutput("plot1.p1", height="320px"), type=6, color = "#E41A1C" ,size= 0.5, hide.ui = FALSE)),
                         column(4, tabsetPanel(
                             tabPanel("Result", 
                                      br(),
                                      span(textOutput("H0.1"), style="color:#318481"), #font-style:italic"
                                      br(),
                                      textOutput("cv1.p1"),
                                      textOutput("pval1.p1")),
                             tabPanel("Summary", 
                                      uiOutput("summary1.p1")
                                      ))
                         )),
                     fluidRow(
                         column(8, 
                                withSpinner(plotOutput("plot2.p1", height="320px"), type = 6, size= 0.5, color = "#E41A1C", hide.ui = FALSE)),
                         column(4, tabsetPanel(
                             tabPanel("Result", 
                                      br(),
                                      span(textOutput("H0.2"), style="color:#318481"),
                                      br(),
                                      textOutput("cv2.p1"),
                                      textOutput("pval2.p1")),
                             tabPanel("Summary", uiOutput("summary2.p1")))
                         )),
                     fluidRow(
                         column(8, withSpinner(plotOutput("plot3.p1", height="320px"), type = 6, size= 0.5, color = "#E41A1C", hide.ui = FALSE)),
                         column(4, tabsetPanel(
                             tabPanel("Result", 
                                      br(),
                                      span(textOutput("H0.3"), style="color:#318481"),
                                      br(),
                                      textOutput("cv3.p1"),
                                      textOutput("pval3.p1")),
                             tabPanel("Summary", uiOutput("summary3.p1")),
                             tabPanel("Interaction Test", tableOutput("intertest.p1")))
                         ))
                     ),
            tabPanel(title="P2 Procedure", 
                     fluidRow(
                         column(8, withSpinner(plotOutput("plot1.p2", height="320px"), type = 6, size= 0.5, color = "#E41A1C", hide.ui = FALSE)),
                         column(4, tabsetPanel(
                             tabPanel("Result", 
                                      br(), br(), br(),
                                      textOutput("cv1.p2"),
                                      textOutput("pval1.p2")),
                             tabPanel("Summary", uiOutput("summary1.p2")))
                         )),
                     fluidRow(
                         column(8, withSpinner(plotOutput("plot2.p2", height="320px"), type = 6, size= 0.5, color = "#E41A1C", hide.ui = FALSE)),
                         column(4, tabsetPanel(
                             tabPanel("Result", 
                                      br(), br(), br(),
                                      textOutput("cv2.p2"),
                                      textOutput("pval2.p2")),
                             tabPanel("Summary", uiOutput("summary2.p2")))
                         )),
                     fluidRow(
                         column(8, withSpinner(plotOutput("plot3.p2", height="320px"), type = 6, size= 0.5, color = "#E41A1C", hide.ui = FALSE)),
                         column(4, tabsetPanel(
                             tabPanel("Result", 
                                      br(), br(), br(),
                                      textOutput("cv3.p2"),
                                      textOutput("pval3.p2")),
                             tabPanel("Summary", uiOutput("summary3.p2")),
                             tabPanel("Interaction Test", tableOutput("intertest.p2")))
                         ))
            )
        )
    )
)



# Server ------------------------------------------------------------------

server <- function(input, output, session) {
    
    output$sampledata<- renderTable( data.frame(`Dependent Variable`= c("Value1", "Value2", "Value3","..."), `Factor A`= c("A1", "A2", "A3", "..."), `Factor B`= c("B1","B2", "B3", "..."), check.names = FALSE) )
    time<- 4
    
    # Read Data
    dt<- eventReactive(input$run,{

        # upload data
        req(input$upload) 
        data<- fread(input$upload$datapath) %>% as.data.frame()
        
        data[,1]<- data[,1] %>% as.numeric
        data[,2]<- data[,2] %>% as.factor
        data[,3]<- data[,3] %>% as.factor
        newdata<- data.frame(data, inter= paste(data[,2],data[,3], sep=":"))  #Define interaction column
        cn<- colnames(data)
        colnames(newdata)[4]<- paste(cn[2], cn[3], sep="*")
        
        return(newdata)
    })
    
    # First Dependency: if data is uploaded, run simulation 
    sim<- eventReactive(dt(), {
      
      showModal(modalDialog("Running in progress", footer=NULL))
      
      # Calculate Statistics
      p1.mean<- lapply(list(dt()[,c(1,2)], dt()[,c(1,3)], dt()[,c(1,4)]), MeanP1)
      p2.mean<- lapply(list(dt()[,c(1,2)], dt()[,c(1,3)], dt()[,c(1,4)]), MeanP2)
  
      # Simulate distribution
      iter<- 10000
      ni <- lapply(p1.mean, function(x) x['ni',])
      n0 <- lapply(p2.mean, function(x) x['n0',])
      distP1<- lapply(1:3, function(i) sapply(1:time, function(q) hdist(ni[[i]], iter)))
      
      removeModal()
      showModal(modalDialog("Still running...", footer=NULL))
      
      distP2<- lapply(1:3, function(i) sapply(1:time, function(q) ddist(n0[[i]], iter)))
      
      removeModal()
      return(list(p1.mean= p1.mean, p2.mean= p2.mean, distP1= distP1, distP2= distP2))
  
    })
    
    
    # Second Dependency: if simulation result is found, calculate critical value and p-value
    test<- eventReactive(c(input$alpha,sim()), {
      #showModal(modalDialog("Updating details", footer=NULL))
      
      cv.p1<-sapply(1:3, function(i) sapply(1:time, function(q) quantile(sim()$distP1[[i]][2,q]$Max, 1-input$alpha/2)) )
      pval.p1<- sapply(1:3, function(i) sapply(1:time, function(q) pvalP1(sim()$p1.mean[[i]], sim()$distP1[[i]][2,q]$Max)))
      out.p1<- cbind(cvmean= apply(cv.p1,2, mean), cvse = apply(cv.p1,2, sd),pval = apply(pval.p1,2, mean))
      
      #removeModal()
      #showModal(modalDialog("Almost there...", footer=NULL))
      
      cv.p2<-sapply(1:3, function(i) apply(sim()$distP2[[i]], 2, function(q) quantile(q, 1-input$alpha/2)) )
      pval.p2<- sapply(1:3, function(i) sapply(1:time, function(q) pvalP2(sim()$p2.mean[[i]], sim()$distP2[[i]][,q])))
      out.p2<- cbind(cvmean= apply(cv.p2,2, mean), cvse = apply(cv.p2,2, sd),pval = apply(pval.p2,2, mean))
      
      #removeModal()
      
      return(list(outP1= out.p1, outP2= out.p2))
      
    })
    
    # Third Dependency: if cv and p-value is calculated, update statistics and charts
    
    out.chartval<- eventReactive(test(),{
      showModal(modalDialog("Cool stuff coming up...", footer=NULL))
      
      v1.p1<- MetricsP1(sim()$p1.mean[[1]], test()$outP1[1,1]) 
      v2.p1<- MetricsP1(sim()$p1.mean[[2]], test()$outP1[2,1] )
      v3.p1<- MetricsP1(sim()$p1.mean[[3]], test()$outP1[3,1])
      
      v1.p2<- MetricsP2(sim()$p2.mean[[1]], test()$outP2[1,1]) 
      v2.p2<- MetricsP2(sim()$p2.mean[[2]], test()$outP2[2,1] )
      v3.p2<- MetricsP2(sim()$p2.mean[[3]], test()$outP2[3,1])
      
      removeModal()
      
      return(list(v1.p1=v1.p1, v2.p1=v2.p1, v3.p1=v3.p1, v1.p2=v1.p2, v2.p2= v2.p2, v3.p2=v3.p2))
    })
    
    # Generate Plots 
    varname<- reactive(colnames(dt()))
    output$plot1.p1<- renderPlot(ggChart2(out.chartval()$v1.p1, text1= paste0("Main Effect: ",varname()[2]), text2= "Treatment mean outside of decision limits suggests statistical difference from the overall average.", text3= "" ))
    output$plot2.p1<- renderPlot(ggChart2(out.chartval()$v2.p1, text1= paste0("Main Effect: ",varname()[3]), text2= "Treatment mean outside of decision limits suggests statistical difference from the overall average.", text3= "" ))
    output$plot3.p1<- renderPlot(ggChart2(out.chartval()$v3.p1, text1= paste0("Interaction Effect: ",varname()[4]), text2= "Treatment mean outside of decision limits suggests statistical difference from the overall average.", text3= "" ) + 
                                   labs(caption = "The colon operator ':' indicates interaction among two or more treatments")
    )
    
    output$plot1.p2<- renderPlot(ggChart2(out.chartval()$v1.p2, text1= paste0("Main Effect: ",varname()[2]), text2= "Treatment mean outside of decision limits suggests statistical difference from the overall average.", text3= "" ))
    output$plot2.p2<- renderPlot(ggChart2(out.chartval()$v2.p2, text1= paste0("Main Effect: ",varname()[3]), text2= "Treatment mean outside of decision limits suggests statistical difference from the overall average.", text3= "" ))
    output$plot3.p2<- renderPlot(ggChart2(out.chartval()$v3.p2, text1= paste0("Interaction Effect: ",varname()[4]), text2= " ", text3= "" ) + 
                                   labs(caption = "The colon operator ':' indicates interaction among two or more treatments"))
    
    
    # Print Summary Table
    output$v1p1 <- renderTable({
      out<- out.chartval()$v1.p1
      rownames(out) <- c("\\( n_{i} \\)", "\\( \\bar{X_i} \\)", "\\( S_{i} \\)", "\\( U_{i} \\)", "\\( V_{i} \\)", "\\( \\tilde{X_i} \\)", "Center", "LDL", "UDL")
      out
      },include.rownames = TRUE, include.colnames = TRUE, striped=TRUE, hover= TRUE, align= "c")
    output$v2p1 <- renderTable({
      out<- out.chartval()$v2.p1
      rownames(out) <- c("\\( n_{i} \\)", "\\( \\bar{X_i} \\)", "\\( S_{i} \\)", "\\( U_{i} \\)", "\\( V_{i} \\)", "\\( \\tilde{X_i} \\)", "Center", "LDL", "UDL")
      out
    },include.rownames = TRUE, include.colnames = TRUE, striped=TRUE, hover= TRUE, align= "c")
    output$v3p1 <- renderTable({
      out<- out.chartval()$v3.p1
      rownames(out) <- c("\\( n_{i} \\)", "\\( \\bar{X_i} \\)", "\\( S_{i} \\)", "\\( U_{i} \\)", "\\( V_{i} \\)", "\\( \\tilde{X_i} \\)", "Center", "LDL", "UDL")
      out
    },include.rownames = TRUE, include.colnames = TRUE, striped=TRUE, hover= TRUE, align= "c")
   
    output$v1p2 <- renderTable({
      out<- out.chartval()$v1.p2
      rownames(out) <- c("\\( n_{i} \\)", "\\( \\bar{X_i} \\)", "\\( S_{i} \\)", "\\( U_{i} \\)", "\\( V_{i} \\)", "\\( \\tilde{X_i} \\)", "Center", "LDL", "UDL")
      out
    },include.rownames = TRUE, include.colnames = TRUE, striped=TRUE, hover= TRUE, align= "c")
    output$v2p2 <- renderTable({
      out<- out.chartval()$v2.p2
      rownames(out) <- c("\\( n_{i} \\)", "\\( \\bar{X_i} \\)", "\\( S_{i} \\)", "\\( U_{i} \\)", "\\( V_{i} \\)", "\\( \\tilde{X_i} \\)", "Center", "LDL", "UDL")
      out
    },include.rownames = TRUE, include.colnames = TRUE, striped=TRUE, hover= TRUE, align= "c")
    output$v3p2 <- renderTable({
      out<- out.chartval()$v3.p2
      rownames(out) <- c("\\( n_{i} \\)", "\\( \\bar{X_i} \\)", "\\( S_{i} \\)", "\\( U_{i} \\)", "\\( V_{i} \\)", "\\( \\tilde{X_i} \\)", "Center", "LDL", "UDL")
      out
    },include.rownames = TRUE, include.colnames = TRUE, striped=TRUE, hover= TRUE, align= "c")
    
  
    ## LaTeX Mode:
    output$summary1.p1 <- renderUI({
      input$alpha # in order to re-render when input$alpha changes
      tagList(
        withMathJax(),
        withMathJax(tableOutput("v1p1"))
      )
    })
    output$summary2.p1 <- renderUI({
      input$alpha # in order to re-render when input$alpha changes
      tagList(
        withMathJax(),
        withMathJax(tableOutput("v2p1"))
      )
    })
    output$summary3.p1 <- renderUI({
      input$alpha # in order to re-render when input$alpha changes
      tagList(
        withMathJax(),
        withMathJax(tableOutput("v3p1"))
      )
    })
    output$summary1.p2 <- renderUI({
      input$alpha # in order to re-render when input$alpha changes
      tagList(
        withMathJax(),
        withMathJax(tableOutput("v1p2"))
      )
    })
    output$summary2.p2 <- renderUI({
      input$alpha # in order to re-render when input$alpha changes
      tagList(
        withMathJax(),
        withMathJax(tableOutput("v2p2"))
      )
    })
    output$summary3.p2 <- renderUI({
      input$alpha # in order to re-render when input$alpha changes
      tagList(
        withMathJax(),
        withMathJax(tableOutput("v3p3"))
      )
    })


    # Print Results Table
    output$H0.1<- renderText(paste0("Null Hypothesis: All treatment means of main effect (",varname()[2],") are equal"))
    output$H0.2<- renderText(paste0("Null Hypothesis: All treatment means of main effect (",varname()[3],") are equal"))
    output$H0.3<- renderText(paste0("Null Hypothesis: All treatment means of interaction effect (",varname()[4],") are equal"))
    
    output$cv1.p1<- renderText(paste0("Critical value: ", test()$outP1[1,1] %>% sprintf("%.3f",.)," (s.e.: ", test()$outP1[1,2] %>% sprintf("%.3f",.), ")" ))
    output$pval1.p1<- renderText(paste0("p-value: ", test()$outP1[1,3] %>% sprintf("%.3f",.)))
    output$cv2.p1<- renderText(paste0("Critical value: ", test()$outP1[2,1] %>% sprintf("%.3f",.)," (s.e.: ", test()$outP1[2,2] %>% sprintf("%.3f",.), ")" ))
    output$pval2.p1<- renderText(paste0("p-value: ", test()$outP1[2,3] %>% sprintf("%.3f",.)))
    output$cv3.p1<- renderText(paste0("Critical value: ", test()$outP1[3,1] %>% sprintf("%.3f",.)," (s.e.: ", test()$outP1[3,2] %>% sprintf("%.3f",.), ")" ))
    output$pval3.p1<- renderText(paste0("p-value: ", test()$outP1[3,3] %>% sprintf("%.3f",.)))
    
    output$cv1.p2<- renderText(paste0("Critical value: ", test()$outP2[1,1] %>% sprintf("%.3f",.)," (s.e.: ", test()$outP2[1,2] %>% sprintf("%.3f",.), ")" ))
    output$pval1.p2<- renderText(paste0("p-value: ", test()$outP2[1,3] %>% sprintf("%.3f",.)))
    output$cv2.p2<- renderText(paste0("Critical value: ", test()$outP2[2,1] %>% sprintf("%.3f",.)," (s.e.: ", test()$outP2[2,2] %>% sprintf("%.3f",.), ")" ))
    output$pval2.p2<- renderText(paste0("p-value: ", test()$outP2[2,3] %>% sprintf("%.3f",.)))
    output$cv3.p2<- renderText(paste0("Critical value: ", test()$outP2[3,1] %>% sprintf("%.3f",.)," (s.e.: ", test()$outP2[3,2] %>% sprintf("%.3f",.), ")" ))
    output$pval3.p2<- renderText(paste0("p-value: ", test()$outP2[3,3] %>% sprintf("%.3f",.)))
    
    output$intertest.p1<- renderTable(interaction_test(dt()[,1:3])$result, rownames=TRUE, align= "c" )
    output$intertest.p2<- renderTable(interaction_test(dt()[,1:3])$result, rownames=TRUE, align= "c" )
    
}




# Run the application 
shinyApp(ui = ui, server = server)

