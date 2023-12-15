library(shiny)
library(sjPlot)
library(corrplot)
library(dplyr)
library(vtable)
source("./script/functions.R")

ui<- fluidPage(
navbarPage("User Interface:",
           tabPanel("Upload",titlePanel("Uploading Files"),
           sidebarLayout(
           sidebarPanel(
           fileInput("file1", "Choose CSV File"),
           
           # Input: Checkbox if file has header ----
           checkboxInput("header", "Header", TRUE),
           
           # Input: Select separator ----
           radioButtons("sep", "Separator",choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),
                        selected = ","),
           
           strong("Note"),
           helpText("Outputs will be based on,", 
                    "1=cateogry +;",
                    "0=cateory -"),),
           
           mainPanel(
                     strong("Dataset overview"),
                     tableOutput("contents"),
                     br(),
                     strong("Marginal probability distribution for the observers (counts)"),
                     verbatimTextOutput("summary"),
                     strong("Marginal probability distribution for the observers (proportions)"),
                     verbatimTextOutput("summary2")
                                          ))), 
           
           tabPanel("Several observers",
                    titlePanel("Agreement between several observers"),
                    sidebarLayout(
                      sidebarPanel(
                        checkboxGroupInput("variables", label = ""),
                        hr(),
                        radioButtons("stat_many", strong("Select the agreement index"),
                                   choices = list("Proportion of agreement" = 1, "Positive agreement" = 2,
                                                   "Jaccard index" = 3,"Intraclass kappa" = 4,
                                                  "Cohen's kappa" = 5),selected = 1)),
                     mainPanel(

                      conditionalPanel(
                        condition = "input.stat_many == 1",
                        strong("Proportion of agreement"),
                        htmlOutput("sum_pom"),
                        br(),
                        strong("Proportion of agreement between pairs of observers"),
                        plotOutput("plotpair1",width = "75%")),
                      conditionalPanel(
                        condition = "input.stat_many == 2",
                        strong("Positive agreement"),
                        htmlOutput("sum_pam"),
                        br(),
                        strong("Positive agreement between pairs of observers"),
                        plotOutput("plotpair2",width = "75%")),
                      conditionalPanel(
                        condition = "input.stat_many == 3",
                        strong("Jaccard index"),
                        htmlOutput("sum_pcam"),
                        br(),
                        strong("Jaccard index between pairs of observers"),
                        plotOutput("plotpair3",width = "75%")),
                      conditionalPanel(
                        condition = "input.stat_many == 4",
                        strong("Intraclass kappa"),
                        htmlOutput("sum_ki"),
                        br(),
                        strong("Intraclass (Fleiss) kappa between pairs of observers"),
                        plotOutput("plotpair4"),width="75%"),
                      conditionalPanel(
                        condition = "input.stat_many == 5",
                        strong("Cohen's kappa"),
                        htmlOutput("sum_kc"),
                        br(),
                        strong("Kappa (Hubert, Conger) between pairs"),
                        plotOutput("plotpair5",width="75%"))
                    ))),

           tabPanel("Two observers",
                      titlePanel("Agreement between two observers"),
                      sidebarLayout(
                       sidebarPanel(
                        selectInput("row", "Select rater 1", character(0)),
                        selectInput("col", "Select rater 2", character(0)),
                        hr(),
                        radioButtons("stat_two", strong("Select the agreement index"),
                                     choices = list("Proportion of agreement" = 1, "Positive agreement" = 2,
                                                    "Jaccard index" = 3,"Intraclass kappa" = 4,
                                                    "Cohen's kappa" = 5),selected = 1)),
                         mainPanel(
                           strong("Contingency table (with cell %)"),
                           htmlOutput("contingency"),
                           br(),
                           strong("Agreement plot"),
                           plotOutput("plotagree",width = "75%"),
                           br(),
                           conditionalPanel(
                             condition = "input.stat_two == 1",
                             strong("Proportion of agreement"),
                             htmlOutput("sum_po2")),
                           conditionalPanel(
                             condition = "input.stat_two == 2",
                             strong("Positive agreement"),
                             htmlOutput("sum_pa2")),
                           conditionalPanel(
                             condition = "input.stat_two == 3",
                             strong("Jaccard index"),
                             htmlOutput("sum_pca2")),
                           conditionalPanel(
                             condition = "input.stat_two==4",
                             strong("Intraclass kappa"),
                             htmlOutput("sum_ki2")),
                           conditionalPanel(
                             condition = "input.stat_two==5",
                             strong("Cohen's kappa"),
                             htmlOutput("sum_kc2"))
                           ))),
           # # tabPanel("Bootstrap",
           # #          titlePanel("Agreement between several observers"),
           # #          sidebarLayout(
           # #            sidebarPanel(
           # #              checkboxGroupInput("var.boot", label = ""), 
           # #              hr(),
           # #              radioButtons("stat_boot", strong("Select the agreement index"),
           # #                           choices = list("Proportion of agreement" = 1, "Positive agreement" = 2,
           # #                                          "Chamberlain positive agreement" = 3,"Intraclass kappa" = 4,
           # #                                          "Cohen's kappa" = 5),selected = 1)),
           #            # mainPanel(
           #            #   
           #            #   conditionalPanel(
           #            #     condition = "input.stat_boot == 1",
           #            #     strong("Proportion of agreement"),
           #            #     htmlOutput("boot_pom"),
           #            #     br(),
           #            #     strong("Proportion of agreement between pairs of observers"), 
           #            #     plotOutput("plotsta1",width = "75%")),
           #            #   conditionalPanel(
           #            #     condition = "input.stat_boot == 2",
           #            #     strong("Positive agreement"),
           #            #     htmlOutput("boot_pam"),
           #            #     br(),
           #            #     strong("Positive agreement between pairs of observers"),
           #            #     plotOutput("plotsta2",width = "75%")),
           #            #   conditionalPanel(
           #            #     condition = "input.stat_boot == 3",
           #            #     strong("Chamberlain positive agreement"),
           #            #     htmlOutput("boot_pcam"),
           #            #     br(),
           #            #     strong("Chamberlain positive agreement between pairs of observers"),
           #            #     plotOutput("plotsta3",width = "75%")),
           #            #   conditionalPanel(
           #            #     condition = "input.stat_boot == 4",
           #            #     strong("Intraclass kappa"),
           #            #     htmlOutput("boot_ki"),
           #            #     br(),
           #            #     strong("Intraclass (Fleiss) kappa between pairs of observers"),
           #            #     plotOutput("plotsta4"),width="75%"),
           #            #   conditionalPanel(
           #            #     condition = "input.stat_boot == 5",
           #            #     strong("Kappa"),
           #            #     htmlOutput("boot_kc"),
           #            #     br(),
           #            #     strong("Kappa between pairs"),
           #            #     plotOutput("plotsta5",width="75%"))
           #            # ))) ,
           tabPanel("Bayesian analysis",
                    titlePanel("Bayesian analysis (Perk's priors)"),
                    sidebarLayout(
                      sidebarPanel(
                        checkboxGroupInput("var_bayes", label = ""),
                        hr(),
                        radioButtons("stat_bayes", strong("Select the agreement index"),
                                     choices = list("Proportion of agreement" = 1, "Positive agreement" = 2,
                                                    "Jaccard index" = 3,"Intraclass kappa" = 4,"Cohen kappa" = 5),selected = 1)),
                      mainPanel(
                        h4("Posterior distribution"),
                        conditionalPanel(
                          condition = "input.stat_bayes == 1", htmlOutput("bayes_pom")),
                        conditionalPanel(
                          condition = "input.stat_bayes == 2",htmlOutput("bayes_pam")),
                        conditionalPanel(
                          condition = "input.stat_bayes== 3", htmlOutput("bayes_pcam")),
                         conditionalPanel(
                          condition = "input.stat_bayes== 4", htmlOutput("bayes_kim")),
                         conditionalPanel(
                          condition = "input.stat_bayes== 5", htmlOutput("bayes_kcm"))
                      ))),

           tabPanel("Sample size",
                    titlePanel("Simple size calculation for varying number of raters"),
                    sidebarLayout(
                      sidebarPanel(
                        radioButtons("stat_size", strong("Select the agreement index"),
                                     choices = list("Proportion of agreement" = 1, "Positive agreement" = 2,
                                                    "Jaccard index" = 3,"Intraclass kappa (Rotondi)" = 4,"Intraclass kappa (Delta)" = 5
                        ),selected = 4),
                        sliderInput("nrater", label = "Number of raters:",
                          min = 2, value = c(3,10), max = 30),
                        sliderInput("nprev", label = "target P(+):",
                                    min = 0.05, value = 0.5, max = 0.95),
                        sliderInput("nagree", label = "target agreement:",
                                    min = 0.1, value = 0.8, max = 1),
                        sliderInput("nwidth", label = "target width:",
                                    min = 0.05, value = 0.1, max = 0.3),
                        sliderInput("ncover", label = "Scenario covered (%):",
                                    min = 50, value = 100, max = 100),
                        sliderInput("ntolo", label = "tolerance agreement:",
                                   min = 0, value = 0, max = 0.1),
                        sliderInput("ntole", label = "tolerance P(+):",
                                   min = 0, value = 0, max = 0.1)
                        ),
                    mainPanel(
                              conditionalPanel(
                                condition="input.stat_size==1",
                                strong("Proportion of agreement"),
                                htmlOutput("n_po")),
                              conditionalPanel(
                                condition ="input.stat_size==2",
                                strong("Positive agreement"),
                                htmlOutput("n_pa")),
                              conditionalPanel(
                                condition = "input.stat_size==3",
                                strong("Jaccard index"),
                                htmlOutput("n_pca")),
                              conditionalPanel(
                                condition= "input.stat_size==4",
                                strong("Intraclass kappa (Rotondi)"),
                                htmlOutput("n_kappai")),
                              conditionalPanel(
                                condition= "input.stat_size==5",
                                strong("Intraclass kappa (Delta)"),
                                htmlOutput("n_kid"))
                              )
                      )),
           tabPanel("Explanations",
                    titlePanel("Explanations"),

                    mainPanel(
                      h4("Agreement plot"),
                      p("The agreement plot is a visual summary of the cross-classification table.
                        The graph is divided in 4 rectangles according to the marginal probability of the two observers.
                        For example, if observer 1 rates 60% of the participants as positive, an horizontal line will be drawn at 60.
                        For observer 2, a vertical line is drawn. The smallest side of the upper left rectangle represents
                        the maximum proportion of participants on which the two obervers can agree for the category 0.
                        Similarly, the lower right rectangle is relative to category 1."),
                      p("The black rectangles represent the current proportion of agreement on category 0 (upper left) and category 1 (lower right).
                         The maximum possible agreement is obtained when one side of a black rectangle coincides with one side of the white rectangles."),
                      p("When the white rectangles cross on the red line, this means that the two observers rated the same proportion of
                        participants in category 1 (i.e., they have the same marginal probability distribution). If the rectangles
                        cross below the red line, this means that a smaller proportion of participants was rated positive by the first observer
                        than by second observer. The inverse is true when the 4 rectangles cross above the red lines."),
                      br(),
                      h4("Correlogram"),
                      p("The correlogram gives the proportion of agreement (po), of positive agreement (pa) or Chamberlain positive agreement (pca)
                         between all possible pairs of observers. The size of a circle is proportional to the value of the agreement, given
                         in white in the circles. Different colors are also used depending on the agreement level."),
                      p("This plot permits to determine if one particular observer contributes very poorly to the overall agreement index
                         or if the agreement is homogeneous between all pairs of observers."),
                      br(),
                      h4("Agreement between several observers"),
                      p("The agreement indexes are weighted means of the pairwise indexes. Pairs with a higher proportion of participants
                        rated positive have a higher weight in the multi-observer agreement index. They are not the mean of the pairwise agreement indexes,
                        except when all pairs have the same propensity to classify participants as positive. In that case, all weights will be equal.")

           ))
  ))



server <- function(input, output, session) {  
  
  output$po_m<-renderText(input$po_many)

  data <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    df1<-read.csv(inFile$datapath,header = input$header,sep = input$sep)
    updateCheckboxGroupInput(session, "variables", label = "Select the observers", choices = colnames(df1), selected = colnames(df1))
    updateCheckboxGroupInput(session, "var_bayes", label = "Select the observers", choices = colnames(df1), selected = colnames(df1))
    return(df1)
  })
  
  data2 <- reactive({
    req(data(),input$variables)
    df2 <- data() %>% select(input$variables)
    return(df2)
  })
  
  data3 <- reactive({
    req(data(),input$var_bayes)
    df3 <- data() %>% select(input$var_bayes)
    return(df3)
  })
  
  
  observeEvent(data(), {
    updateSelectInput(session, "row", choices = names(data()))
    updateSelectInput(session, "col", choices = names(data()))
  })
    
 
   output$contingency <- renderUI({
      req(data())
      custom_table<-sjt.xtab(data()[[input$row]],data()[[input$col]],var.labels=c(paste(input$row),paste(input$col)),show.cell.prc=TRUE,show.summary=FALSE)
      HTML(custom_table$knitr)
   })

   output$contents <- renderTable({
     req(data())
     head(data())
   })

  output$summary <- renderPrint({
    
    req(data())
     #compute unique levels in data frame
    lvls <- unique(unlist(data()))
    
    # apply the summation per value 
    sapply(data(),function(x) table(factor(x, levels = lvls,ordered = TRUE)))
    
  })
  
  output$summary2 <- renderPrint({
    
    req(data())
    #compute unique levels in data frame
    lvls <- unique(unlist(data()))
    
    # apply the summation per value 
    sapply(data(),function(x) round(prop.table(table(factor(x, levels = lvls,ordered = TRUE))),2))
    
  })
  
  output$plotagree <- renderPlot({

    req(data())
    
    agree.plot(data()[[input$row]],data()[[input$col]])
    
  })
  
  output$plotpair1 <- renderPlot({
    
    req(data2())
    
    correlo1(data2())
    
  })
  
  output$plotpair2 <- renderPlot({
    
    req(data2())
    
    correlo2(data2())
    
  })
  
  output$plotpair3 <- renderPlot({
    
    req(data2())
    
    correlo3(data2())
    
  })
  
  output$plotpair4 <- renderPlot({
    
    req(data2())
    
    correlo4(data2())
    
  })
  
  output$plotpair5 <- renderPlot({
    
    req(data2())
    
    correlo5(data2())
    
  })
  
  
  
  
  
  output$sum_po2 <- renderUI({
    req(data())
    
    data2<-as.data.frame(cbind(data()[[input$row]],data()[[input$col]]))
    
    po<-round(po(data2),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("po"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))
      ))
  })
    
  output$sum_pa2 <- renderUI({
    req(data())
    
    data2<-as.data.frame(cbind(data()[[input$row]],data()[[input$col]]))
    
    pa<-round(pa(data2),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("ps"),
        tags$td(paste(pa[1])),
        tags$td(paste(pa[2])),
        tags$td(paste(paste(paste(round(pa[1]-1.96*pa[2],3)),";")),paste(round(pa[1]+1.96*pa[2],3)))
      ))
  })
    
  
  output$sum_pca2 <- renderUI({
    req(data())
    
    data2<-as.data.frame(cbind(data()[[input$row]],data()[[input$col]]))
    
    pca<-round(pca(data2),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("pJ"),
        tags$td(paste(pca[1])),
        tags$td(paste(pca[2])),
        tags$td(paste(paste(paste(round(pca[1]-1.96*pca[2],3)),";")),paste(round(pca[1]+1.96*pca[2],3)))
      ))
  })
  
  
  output$sum_ki2 <- renderUI({
    req(data())
    
    data2<-as.data.frame(cbind(data()[[input$row]],data()[[input$col]]))
    
    ki<-round(many1(data2),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("KI"),
        tags$td(paste(ki[1])),
        tags$td(paste(ki[2])),
        tags$td(paste(paste(paste(round(ki[1]-1.96*ki[2],3)),";")),paste(round(ki[1]+1.96*ki[2],3)))
      ))
  })
  
  
  output$sum_kc2 <- renderUI({
    req(data())
    
    data2<-as.data.frame(cbind(data()[[input$row]],data()[[input$col]]))
    
    kc<-round(many2(data2),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("Kc"),
        tags$td(paste(kc[1])),
        tags$td(paste(kc[2])),
        tags$td(paste(paste(paste(round(kc[1]-1.96*kc[2],3)),";")),paste(round(kc[1]+1.96*kc[2],3)))
      ))
  })
  
  output$sum_pom <- renderUI({
    req(data2())
    
    po<-round(po(data2()),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("po"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))))
    
  }) 
  
  output$sum_pam <- renderUI({
    req(data2())
    
    pa<-round(pa(data2()),3)
    pa_m<-round(pa_mean(data2()),3)
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval"),
        tags$th("Mean pa")
      ),
      tags$tr(
        tags$td("ps"),
        tags$td(paste(pa[1])),
        tags$td(paste(pa[2])),
        tags$td(paste(paste(paste(round(pa[1]-1.96*pa[2],3)),";")),paste(round(pa[1]+1.96*pa[2],3))),
        tags$td(paste(pa_m))))
    
  }) 
  
   
  
  output$sum_pcam <- renderUI({
    req(data2())
    
    pca<-round(pca(data2()),3)
    pca_m<-round(pca_mean(data2()),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval"),
        tags$th("Mean pca")
      ),
      tags$tr(
        tags$td("pJ"),
        tags$td(paste(pca[1])),
        tags$td(paste(pca[2])),
        tags$td(paste(paste(paste(round(pca[1]-1.96*pca[2],3)),";")),paste(round(pca[1]+1.96*pca[2],3))),
        tags$td(paste(pca_m))))
    
  }) 
  
  output$sum_kc<- renderUI({
    req(data2())
    
    kc<-round(many2(data2()),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("Kc"),
        tags$td(paste(kc[1])),
        tags$td(paste(kc[2])),
        tags$td(paste(paste(paste(round(kc[1]-1.96*kc[2],3)),";")),paste(round(kc[1]+1.96*kc[2],3)))))
    
  }) 
 
  
  
  output$sum_ki<- renderUI({
    req(data2())
    
    ki<-round(many1(data2()),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("KI"),
        tags$td(paste(ki[1])),
        tags$td(paste(ki[2])),
        tags$td(paste(paste(paste(round(ki[1]-1.96*ki[2],3)),";")),paste(round(ki[1]+1.96*ki[2],3)))))
    
  }) 
  
  
  output$n_po <- renderTable({
  
     length<-input$nrater[2]-input$nrater[1]+1
     count<-1
     n<-rep(NA,length)
     for (r in c(input$nrater[1]:input$nrater[2])){
    n[count]<-size.pom(w=input$nwidth,pm=input$nprev,po=input$nagree,a.level=0.05,e_po=input$ntolo,e_pm=input$ntole,quant=input$ncover,nrat=r,nb_s=40000)
    count<-count+1
     }
     res.n<-data.frame(round(cbind(seq(input$nrater[1],input$nrater[2]),n),0))
     return(res.n)
    
  })
  
  
  output$n_pa <- renderTable({
    
    length<-input$nrater[2]-input$nrater[1]+1
    count<-1
    n<-rep(NA,length)
    for (r in c(input$nrater[1]:input$nrater[2])){
      n[count]<-size.pam(w=input$nwidth,pm=input$nprev,pa=input$nagree,a.level=0.05,e_pa=input$ntolo,e_pm=input$ntole,quant=input$ncover,nrat=r,nb_s=40000)
      count<-count+1
    }
    res.n<-data.frame(round(cbind(seq(input$nrater[1],input$nrater[2]),n),0))
    return(res.n)
    
  })
  
  
  output$n_pca <- renderTable({
    
    length<-input$nrater[2]-input$nrater[1]+1
    count<-1
    n<-rep(NA,length)
    for (r in c(input$nrater[1]:input$nrater[2])){
      n[count]<-size.pcam(w=input$nwidth,pm=input$nprev,pca=input$nagree,a.level=0.05,e_pca=input$ntolo,e_pm=input$ntole,quant=input$ncover,nrat=r,nb_s=40000)
      count<-count+1
    }
    res.n<-data.frame(round(cbind(seq(input$nrater[1],input$nrater[2]),n),0))
    return(res.n)
    
  })

  output$n_kappai <- renderTable({
    
    length<-input$nrater[2]-input$nrater[1]+1
    count<-1
    n<-rep(NA,length)
    for (r in c(input$nrater[1]:input$nrater[2])){
      n[count]<-size.kappai(w=input$nwidth,prop=input$nprev,kappa0=input$nagree,raters=r)
      count<-count+1
    }
    res.n<-data.frame(round(cbind(seq(input$nrater[1],input$nrater[2]),n),0))
    colnames(res.n)<-c("nraters","nsubjects")
    return(res.n)
    
  })
  
  output$n_kid <- renderTable({
    
    length<-input$nrater[2]-input$nrater[1]+1
    count<-1
    n<-rep(NA,length)
    for (r in c(input$nrater[1]:input$nrater[2])){
      n[count]<-size.kid(w=input$nwidth,pm=input$nprev,ki=input$nagree,a.level=0.05,e_ki=input$ntolo,e_pm=input$ntole,quant=input$ncover,nrat=r,nb_s=40000)
      count<-count+1
    }
    res.n<-data.frame(round(cbind(seq(input$nrater[1],input$nrater[2]),n),0))
    return(res.n)
    
  })
  # output$boot_pom <- renderUI({
  #   req(data3())
  #   set.seed(input$seed)
  #   po.sample<-boot(data3(), po.boot,input$nboot)
  #   
  #   tags$table(
  #     style = "width:100%",
  #     tags$tr(
  #       tags$th(""),
  #       tags$th("Value"),
  #       tags$th("95% confidence interval")
  #     ),
  #     tags$tr(
  #       tags$td("po"),
  #       tags$td(paste(round(quantile(po.sample, probs = 0.5),3))),
  #       tags$td(paste(paste(paste(round(quantile(po.sample, probs = 0.025),3)),";")),paste(round(quantile(po.sample, probs = 0.975),3)))))
  # }) 
  
  
  output$bayes_pom <- renderUI({
    req(data3())
    
    po_b<-round(po.bayes(data3(),c(0.025,0.50,0.975)),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Mean"),
        tags$th("SE"),
        tags$th("Median"),
        tags$th("95% credibility interval (non-parametric)")
      ),
      tags$tr(
        tags$td("po"),
        tags$td(paste(po_b[1])),
        tags$td(paste(po_b[2])),
        tags$td(paste(po_b[4])),
        tags$td(paste(paste(paste(po_b[3]),";"),paste(po_b[5])))))
    
  }) 
  
  output$bayes_pam <- renderUI({
    req(data3())
    
    pa_b<-round(pa.bayes(data3(),5000,c(0.025,0.50,0.975)),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Mean"),
        tags$th("SE"),
        tags$th("Median"), 
        tags$th("95% credibility interval (non-parametric)")
      ),
      tags$tr(
        tags$td("pa"),
        tags$td(paste(pa_b[1])),
        tags$td(paste(pa_b[2])),
        tags$td(paste(pa_b[4])),    
        tags$td(paste(paste(paste(pa_b[3]),";"),paste(pa_b[5])))))
    
  }) 
  
  output$bayes_pcam <- renderUI({
    req(data3())
    
    pca_b<-round(pca.bayes(data3(),5000,c(0.025,0.50,0.975)),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Mean"),
        tags$th("SE"),
        tags$th("Median"),
        tags$th("95% credibility interval (non-parametric)")
      ),
      tags$tr(
        tags$td("pca"),
        tags$td(paste(pca_b[1])),
        tags$td(paste(pca_b[2])),
        tags$td(paste(pca_b[4])),
        tags$td(paste(paste(paste(pca_b[3]),";"),paste(pca_b[5])))))
    
  }) 
  
  output$bayes_kcm <- renderUI({
    req(data3())
    
    kc_b<-round(kc.bayes(data3(),c(0.025,0.50,0.975)),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Mean"),
        tags$th("SE"),
        tags$th("Median"),
        tags$th("95% credibility interval (non-parametric)")
      ),
      tags$tr(
        tags$td("kc"),
        tags$td(paste(kc_b[1])),
        tags$td(paste(kc_b[2])),
        tags$td(paste(kc_b[4])),
        tags$td(paste(paste(paste(kc_b[3]),";"),paste(kc_b[5])))))
    
  }) 
  
  
  output$bayes_kim<- renderUI({
    req(data3())
    
    ki_b<-round(ki.bayes(data3(),c(0.025,0.50,0.975)),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Mean"),
        tags$th("SE"),
        tags$th("Median"),
        tags$th("95% credibility interval (non-parametric)")
      ),
      tags$tr(
        tags$td("ki"),
        tags$td(paste(ki_b[1])),
        tags$td(paste(ki_b[2])),
        tags$td(paste(ki_b[4])),
        tags$td(paste(paste(paste(ki_b[3]),";"),paste(ki_b[5])))))
    
  }) 
  
}

shinyApp(ui, server)