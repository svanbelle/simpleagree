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
                    "1=positive;",
                    "0=negative"),),
           
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
                                                   "Chamberlain positive agreement" = 3),selected = 1)),
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
                        strong("Chamberlain positive agreement"),
                        htmlOutput("sum_pcam"),
                        br(),
                        strong("Chamberlain positive agreement between pairs of observers"),
                        plotOutput("plotpair3",width = "75%"))
                    ))) ,
           
           tabPanel("Two observers",
                      titlePanel("Agreement between two observers"),
                      sidebarLayout(
                       sidebarPanel(
                        selectInput("row", "Select rater 1", character(0)),
                        selectInput("col", "Select rater 2", character(0)),
                        hr(),
                        radioButtons("stat_two", strong("Select the agreement index"),
                                     choices = list("Proportion of agreement" = 1, "Positive agreement" = 2,
                                                    "Chamberlain positive agreement" = 3),selected = 1)),
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
                             strong("Chamberlain positive agreement"),
                             htmlOutput("sum_pca2"))
                           ))),
           
           tabPanel("Bayesian analysis",
                    titlePanel("Bayesian analysis (Perk's priors)"),
                    sidebarLayout(
                      sidebarPanel(
                        checkboxGroupInput("var_bayes", label = ""), 
                        hr(),
                        radioButtons("stat_bayes", strong("Select the agreement index"),
                                     choices = list("Proportion of agreement" = 1, "Positive agreement" = 2,
                                                    "Chamberlain positive agreement" = 3),selected = 1)),
                    mainPanel(
                    h4("Posterior distribution"),
                    conditionalPanel(
                      condition = "input.stat_bayes == 1", htmlOutput("bayes_pom")),
                    conditionalPanel(
                      condition = "input.stat_bayes == 2",htmlOutput("bayes_pam")),
                    conditionalPanel(
                      condition = "input.stat_bayes== 3", htmlOutput("bayes_pcam")) 
                      
                    ))),
           
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
        tags$td("pa"),
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
        tags$td("pca"),
        tags$td(paste(pca[1])),
        tags$td(paste(pca[2])),
        tags$td(paste(paste(paste(round(pca[1]-1.96*pca[2],3)),";")),paste(round(pca[1]+1.96*pca[2],3)))
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
        tags$td("pa"),
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
        tags$td("pca"),
        tags$td(paste(pca[1])),
        tags$td(paste(pca[2])),
        tags$td(paste(paste(paste(round(pca[1]-1.96*pca[2],3)),";")),paste(round(pca[1]+1.96*pca[2],3))),
        tags$td(paste(pca_m))))
    
  }) 
  
  
 
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
  
}

shinyApp(ui, server)