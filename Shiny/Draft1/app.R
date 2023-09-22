#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


pacman::p_load(shiny, data.table, dplyr, shinythemes, waiter, ggplot2, ggtext, ggthemr, waiter, remotes)

# Functions ---------------------------------------------------------------
waiting_screen <- tagList(
    spin_pong(),
    h4("Cool stuff loading...")
) 

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
    Ui <- sapply(1:nk, function(a) (1/ni[a]) + (1/ni[a]) * sqrt((1/(ni[a]-1))*(maxvar/vari[a] - 1)) ) #Eq4
    Vi <- sapply(1:nk, function(a) (1/ni[a]) - (1/ni[a]) * sqrt((ni[a]-1)*(maxvar/vari[a] - 1)) )  #Eq5
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


# Simulate h* distribution (critical values)
OneRunP1<- function(ni){
    k<- length(ni)
    ti<- sapply(1:k, function(p) rt(1, ni[p]-2))   #t(ni-2) distribution
    ttildei<- sapply(1:k, function(p) ((k-1)/k)*ti[p]-(sqrt(ni[p])/k)*sum(ti[-p]/sqrt(ni[-p])))  #Eq9 (transform student-t dis into ttilde)
    cv<- c(Min= min(ttildei), Max= max(ttildei))
    return(cv)
}
OneRunP2<- function(W){
    n0<- W['n0',] %>% unique()
    k<- ncol(W)
    t <- rt(n = k, df = n0-1)
    max(t - mean(t))
}

# Calculate critical values at level alpha, based on cvdist (must be an object from SimCV() ).
CalcCV<- function(cvdist, alpha){
    CV.k<- quantile(cvdist[,2], 1-alpha/2)
    CV.1<- -quantile(cvdist[,1], alpha/2)
    max(CV.1, CV.k)  #Eq11
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
    pquantile<- sapply(pseq, function(w) quantile(cvdist[,2], 1-w/2))
    pequation<- abs(pquantile- max(abs(tval)))
    pval<- pseq[pequation== min(pequation)]
    
    return(pval)
    
    #if(P.value1 > 0 && P.value1 < 1){
    #  PV1<-sprintf("%.3f",P.value1)
    #}else if(P.value1 == 0){
    #  PV1<-"< 0.001 "
    #}else{PV1 <-"> 0.999"}
    
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

# Compute CV and P-value
testP1<- function(wlist, dist, alpha){
    k<- length(dist)
    time<- 4
    out<- NULL
    
    for(q in 1:k){
        h<- sapply(1:time, function(i) CalcCV(dist[[q]][,c(1:2)*i], alpha))
        h<- data.frame(`CV Mean`= mean(h), `CV S.E.`= sd(h))
        
        p.val<- sapply(1:time, function(i) pvalP1(wlist[[q]], dist[[q]][,c(1:2)*i]) ) %>% mean
        h<- cbind(h, pval= p.val) 
        out<- rbind(out, h)
        
    }
    rownames(out) <- c(1:k)
    return(out)
}

testP2<- function(wlist, dist, alpha){
    k<- length(dist)
    time<- 4
    out<- NULL
    
    for(q in 1:k){
        # Calculate critical values
        d<- sapply(1:time, function(i) quantile(dist[[q]][,i], 1-alpha/2))
        h<- data.frame(`CV Mean`= mean(d), `CV S.E.`= sd(d) )
        
        # Calculate p-value
        p.val<- sapply(1:time, function(i) pvalP2(wlist[[q]], dist[[q]][,i]) ) %>% mean
        h<- cbind(h, pval= p.val) 
        out<- rbind(out, h)
    }
    rownames(out) <- c(1:k)
    return(out)
}

MetricsP1<- function(W, criticalvalue){
    xtildei <- W["xtildei",]
    k<- ncol(W)
    ni <- W["ni",]
    z1<- max(W["vari",]) 
    
    center <- mean(xtildei)
    LDL <- mean(xtildei) - criticalvalue * sqrt(z1) / sqrt(ni)
    UDL <- mean(xtildei) + criticalvalue * sqrt(z1)/ sqrt(ni)
    
    return(rbind(W, center=center, LDL=LDL, UDL=UDL))
}

MetricsP2<- function(W, criticalvalue){
    xtildei <- W["xtildei",]
    k<- ncol(W)
    ni <- W["ni",]
    z1<- max(W["vari",]/ni)

    center <- mean(xtildei)
    LDL <- mean(xtildei) - criticalvalue * sqrt(z1)
    UDL <- mean(xtildei) + criticalvalue * sqrt(z1)
    
    return(rbind(W, center=center, LDL=LDL, UDL=UDL))
    
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
        labs(title=text1, subtitle= text2, x= text3, y= "Xtilde") +
        theme(plot.margin= margin(20,50,25,25))
    
    
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


# Example usage:
# ChartP1_ggplot(chartvalue, criticalvalue)

# Example usage:
# ChartP1_ggplot(chartvalue, criticalvalue)


# UI ----------------------------------------------------------------------


ui <- fluidPage(
    #autoWaiter(html= spin_pulsar(), color="black"),
    useWaiter(),
    #waiterShowOnLoad(html=spin_square_circle(), color="black"),
    waiterOnBusy(html= spin_dots()),
    #useHostess(),
    theme = shinytheme("united"),

    # Application title
    titlePanel("Two-way HANOM"),
    helpText("Some text here"),
    
    # Sidebar with a slider input for number of bins 
    sidebarPanel(
        #helpText("A 3-column input: first column must be numerical"),
        fileInput("upload","Upload data", multiple = FALSE),
        useWaitress(),
        actionButton("run","Run Analysis"),
        width=3
    ),
    mainPanel(
        numericInput("alpha","Type I error rate: ", value= 0.05, min=0, max=1, step= 0.01),
        tabsetPanel(
            id="tabs",
            tabPanel(title="P1 Procedure", 
                     fluidRow(
                         column(8, plotOutput("plot1.p1", height="320px")),
                         column(4, tabsetPanel(
                             tabPanel("Result", 
                                      br(), br(), br(),
                                      textOutput("cv1.p1"),
                                      textOutput("pval1.p1")),
                             tabPanel("Summary", tableOutput("summary1.p1")))
                         )),
                     fluidRow(
                         column(8, 
                                plotOutput("plot2.p1", height="320px")),
                         column(4, tabsetPanel(
                             tabPanel("Result", 
                                      br(), br(), br(),
                                      textOutput("cv2.p1"),
                                      textOutput("pval2.p1")),
                             tabPanel("Summary", tableOutput("summary2.p1")))
                         )),
                     fluidRow(
                         column(8, plotOutput("plot3.p1", height="320px")),
                         column(4, tabsetPanel(
                             tabPanel("Result", 
                                      br(), br(), br(),
                                      textOutput("cv3.p1"),
                                      textOutput("pval3.p1")),
                             tabPanel("Summary", tableOutput("summary3.p1")),
                             tabPanel("Interaction Test", tableOutput("intertest.p1")))
                         ))
                     ),
            tabPanel(title="P2 Procedure", 
                     fluidRow(
                         column(8, plotOutput("plot1.p2", height="320px")),
                         column(4, tabsetPanel(
                             tabPanel("Result", 
                                      br(), br(), br(),
                                      textOutput("cv1.p2"),
                                      textOutput("pval1.p2")),
                             tabPanel("Summary", tableOutput("summary1.p2")))
                         )),
                     fluidRow(
                         column(8, plotOutput("plot2.p2", height="320px")),
                         column(4, tabsetPanel(
                             tabPanel("Result", 
                                      br(), br(), br(),
                                      textOutput("cv2.p2"),
                                      textOutput("pval2.p2")),
                             tabPanel("Summary", tableOutput("summary2.p2")))
                         )),
                     fluidRow(
                         column(8, plotOutput("plot3.p2", height="320px")),
                         column(4, tabsetPanel(
                             tabPanel("Result", 
                                      br(), br(), br(),
                                      textOutput("cv3.p2"),
                                      textOutput("pval3.p2")),
                             tabPanel("Summary", tableOutput("summary3.p2")),
                             tabPanel("Interaction Test", tableOutput("intertest.p2")))
                         ))
            )
        )
    )
)



# Server ------------------------------------------------------------------

server <- function(input, output, session) {
    
    #w<- Waiter$new(c("plot1.p1", "plot2.p1", "plot3.p1"))
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
        colnames(newdata)[4]<- paste(cn[2], cn[3], sep=":")
        
        return(newdata)
    })
    

    
    #observeEvent(input$run, {
    #    waiter_show(html = waiting_screen, color = "black")
    #    Sys.sleep(5)
    #    waiter_hide()
    #})
# P1 Procedure ------------------------------------------------------------

    # Calculate Statistics
    p1.mean<- eventReactive(dt(), {
        summary<- lapply(list(dt()[,c(1,2)], dt()[,c(1,3)], dt()[,c(1,4)]), MeanP1)
        return(summary)
    })
    
    #Calculate h distribution
    #p1.hdist<-eventReactive(dt(),{
    #    iter<- 1000
    #    times<- 4
    #    totalloop<- iter*times*3
    #        withProgressWaitress({
    #            htemp<- array(0, dim=c(iter,2*times))
    #            out.hdist<- list(htemp, htemp, htemp)
    #            
    #            for(h in 1:3){
    #                for(j in 1:times){
    #                    for(i in 1:iter){
    #                        out.hdist[[h]][i,(2*j-1):(2*j)] <- OneRunP1(p1.mean()[[h]]['ni',]) #Min, Max
    #                        incProgressWaitress(1/(iter*times*3))
    #                    }
    #                }
    #            }
    #        }, selector = "#run", theme = "overlay-percent", max= totalloop, infinite=FALSE)
    #    return(out.hdist)
    # })
    #
    p1.hdist<-eventReactive(dt(),{
        iter<- 1000
        times<- 4
        totalloop<- iter*times*3
            htemp<- array(0, dim=c(iter,2*times))
            out.hdist<- list(htemp, htemp, htemp)
            
            for(h in 1:3){
                for(j in 1:times){
                    for(i in 1:iter){
                        out.hdist[[h]][i,(2*j-1):(2*j)] <- OneRunP1(p1.mean()[[h]]['ni',]) #Min, Max
                    }
                }
            }
        return(out.hdist)
    })
        
    #Calculate critical value and p-value
    p1.test<- reactive(testP1(p1.mean(), p1.hdist(), input$alpha))
    
    output$cv1.p1<- renderText(paste0("Critical value: ", p1.test()[1,1] %>% round(.,3)," (s.e.: ", p1.test()[1,2] %>% round(.,3), ")" ))
    output$pval1.p1<- renderText(paste0("p-value: ", p1.test()[1,3] %>% round(.,3)))
    output$cv2.p1<- renderText(paste0("Critical value: ", p1.test()[2,1] %>% round(.,3)," (s.e.: ", p1.test()[2,2] %>% round(.,3), ")" ))
    output$pval2.p1<- renderText(paste0("p-value: ", p1.test()[2,3] %>% round(.,3)))
    output$cv3.p1<- renderText(paste0("Critical value: ", p1.test()[3,1] %>% round(.,3)," (s.e.: ", p1.test()[3,2] %>% round(.,3), ")" ))
    output$pval3.p1<- renderText(paste0("p-value: ", p1.test()[3,3] %>% round(.,3)))

    output$intertest.p1<- renderTable(interaction_test(dt()[,c(1:3)])$result, rownames=TRUE, align= "c" )
    
    # Calculate Statistics
    out.chartval<- eventReactive(p1.test(),{
        v1<- MetricsP1(p1.mean()[[1]], p1.test()[1,1]) 
        v2<- MetricsP1(p1.mean()[[2]], p1.test()[2,1] )
        v3<- MetricsP1(p1.mean()[[3]], p1.test()[3,1])
        return(list(v1=v1, v2=v2, v3=v3))
    })
    
    output$summary1.p1 <- renderTable(out.chartval()$v1, rownames=TRUE, striped=TRUE, hover= TRUE, align= "c", )
    output$summary2.p1 <- renderTable(out.chartval()$v2, rownames=TRUE, striped=TRUE, hover= TRUE, bordered= FALSE, align= "c")
    output$summary3.p1 <- renderTable(out.chartval()$v3, rownames=TRUE, striped=TRUE, hover= TRUE, bordered= FALSE, align= "c")
   
    # Plots 
    varname<- reactive(colnames(dt()))
    
    output$plot1.p1<- renderPlot({
        input$run
        Sys.sleep(10)
        ggChart2(out.chartval()$v1, text1= paste0("Main Effect: ",varname()[2]), text2= "HANOM P1 Results: ", text3= varname()[2] )
    })

    output$plot2.p1<- renderPlot(ggChart2(out.chartval()$v2, text1= paste0("Main Effect: ",varname()[3]), text2= "HANOM P1 Results: ", text3= varname()[3] ))
    output$plot3.p1<- renderPlot(ggChart2(out.chartval()$v3, text1= paste0("Main Effect: ",varname()[4]), text2= "HANOM P1 Results: ", text3= varname()[4] ))
    
    

# P2 Procedure ------------------------------------------------------------
    
    # Calculate Statistics
    p2.mean<- eventReactive(dt(), {
        summary<- lapply(list(dt()[,c(1,2)], dt()[,c(1,3)], dt()[,c(1,4)]), MeanP2)
        return(summary)
    })
    
    # Calculate d distribution
    p2.ddist<-eventReactive(dt(),{
        iter<- 100
        times<- 4
        totalloop<- iter*times*3
        dtemp<- array(0, dim=c(iter,times))  
        out.ddist<- list(dtemp, dtemp, dtemp)
        
        for(h in 1:3){
            for(j in 1:times){
                for(i in 1:iter){
                    out.ddist[[h]][i,j] <- OneRunP2(p2.mean()[[h]]) #Min, Max
                }
            }
        }
        return(out.ddist)
        message("yes")
    })
    
    # Calculate critical value and p-value
    p2.test<- reactive(testP2(p2.mean(), p2.ddist(), input$alpha))
    
    output$cv1.p2<- renderText(paste0("Critical value: ", p2.test()[1,1] %>% round(.,3)," (s.e.: ", p2.test()[1,2] %>% round(.,3), ")" ))
    output$pval1.p2<- renderText(paste0("p-value: ", p2.test()[1,3] %>% round(.,3)))
    output$cv2.p2<- renderText(paste0("Critical value: ", p2.test()[2,1] %>% round(.,3)," (s.e.: ", p2.test()[2,2] %>% round(.,3), ")" ))
    output$pval2.p2<- renderText(paste0("p-value: ", p2.test()[2,3] %>% round(.,3)))
    output$cv3.p2<- renderText(paste0("Critical value: ", p2.test()[3,1] %>% round(.,3)," (s.e.: ", p2.test()[3,2] %>% round(.,3), ")" ))
    output$pval3.p2<- renderText(paste0("p-value: ", p2.test()[3,3] %>% round(.,3)))
    
    output$intertest.p2<- renderTable(interaction_test(dt()[,1:3])$result, rownames=TRUE, align= "c" )
    
    # Calculate Statistics
    out.chartval.p2<- eventReactive(p2.test(),{
        v1<- MetricsP2(p2.mean()[[1]], p2.test()[1,1]) 
        v2<- MetricsP2(p2.mean()[[2]], p2.test()[2,1] )
        v3<- MetricsP2(p2.mean()[[3]], p2.test()[3,1])
        return(list(v1=v1, v2=v2, v3=v3))
    })
    
    output$summary1.p2 <- renderTable(out.chartval.p2()$v1, rownames=TRUE, striped=TRUE, hover= TRUE, align= "c", )
    output$summary2.p2 <- renderTable(out.chartval.p2()$v2, rownames=TRUE, striped=TRUE, hover= TRUE, bordered= FALSE, align= "c")
    output$summary3.p2 <- renderTable(out.chartval.p2()$v3, rownames=TRUE, striped=TRUE, hover= TRUE, bordered= FALSE, align= "c")
    
    # Plots 
    output$plot1.p2<- renderPlot(ggChart2(out.chartval.p2()$v1, text1= paste0("Main Effect: ",varname()[2]), text2= "HANOM P2 Results: ", text3= varname()[2] ))
    output$plot2.p2<- renderPlot(ggChart2(out.chartval.p2()$v2, text1= paste0("Main Effect: ",varname()[3]), text2= "HANOM P2 Results: ", text3= varname()[3] ))
    output$plot3.p2<- renderPlot(ggChart2(out.chartval.p2()$v3, text1= paste0("Main Effect: ",varname()[4]), text2= "HANOM P2 Results: ", text3= varname()[4] ))

        
}




# Run the application 
shinyApp(ui = ui, server = server)
