rm(list=ls()); gc()
library(shiny); library(ggplot2)

# Define UI for data upload app ----
ui <- fluidPage(
  # App title ----
  titlePanel("Summarizing Results From Geno-Diver"),
  # Sidebar layout with a input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel
    (
        # take file input
        fluidRow(column(12,wellPanel(fileInput('file1', 'Choose One Of Output Files')))), # END wellPanel()
        # Input: Numeric entry for number of obs to view ----
        # numericInput(inputId = "obs",label = "Number of observations to view:",value = 5),
        # "Empty inputs" - they will be updated after the data is uploaded
        selectInput('ycol', 'Plotting Variable by Generation', "", selected = ""),
        sliderInput(inputId = "startgeneration",label = "Plot Generation (Minimum):",min = 0,max = 0,value = NA,step=1),
        sliderInput(inputId = "endgeneration",label = "Plot Generation (Maximum):",min = 0,max = 0,value = NA,step=1)
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      # Output errors messages if given file that doesn't produce summary statistics
      tags$style(type='text/css', '#text1 {font-size:150%; color: red;}'), 
      textOutput("text1"),
      
      # Output: Formatted text for caption ----
      tags$style(type='text/css', '#caption {font-size:150%;}'), 
      textOutput("caption", container = span),
      # Output Plot of variable by Generation #
      plotOutput('MyPlot')
      # Output: HTML table with requested number of observations ----
      #tableOutput("view")
      
    )
  )
)

server <- function(input, output,session) {
  
  values <- reactiveValues(df_data = NULL)
  COLUMNS <- reactiveValues(df_data = NULL)
  fileinputed <- reactiveValues(filename = "")
  
  ### Read in Datafile ###
  observeEvent(input$file1, 
  {
    # Figure out which file it is #
    tempstring <- unlist(strsplit(as.character(paste(input$file1[[1]])),"_"))
    ## Performance File ##
    datafile = "NotFound";
    ## first check to see if performance file ##
    X <- grep("Performance", tempstring)
    if(length(X) > 0) {datafile = "Performance"}
    ## if not performance check for inbreeding summary stats ##
    if(datafile == "NotFound")
    {
      X <- grep("Inbreeding", tempstring)
      if(length(X) > 0) {datafile = "Inbreeding"}
    }
    if(datafile == "NotFound")
    {
      X <- grep("QTL", tempstring)
      if(length(X) > 0) {datafile = "QTL"}
    }
    if(datafile == "NotFound")
    {
      X <- grep("Decay", tempstring)
      if(length(X) > 0) {datafile = "Decay"}
    }
    ## if still not figured out give an error ##
    if(datafile == "NotFound"){output$text1 <- renderText({"No Summary Statistics for DataSet"})}
    ## if figure out read in dataframe and put it into right format ##
    if(datafile != "NotFound")
    {
      # read the file in
      if (is.null(input$file1$datapath)) return(NULL)
      # read in file
      if(datafile == "Performance")
      {
        fileinputed$filename = "Performance";
        values$df_data <- read.table(input$file1$datapath, header=TRUE, sep = " ",colClasses = c("numeric",rep("character",6)))
        # Remove Standard Deviations #
        temp <- values$df_data 
        for(i in 2:ncol(temp)){temp[ ,i] <- as.numeric(matrix(unlist(strsplit(temp[,i],'[()]')),ncol=2,byrow=TRUE)[,1])}
        temp <- apply(temp,2,as.numeric)
        values$df_data <- data.frame(temp,stringsAsFactors = FALSE);
        ## Update columns that can be plotted ##
        variables <- names(values$df_data)
        variables <- variables[2:length(variables)]
        COLUMNS$df_data <- c("Generation","Phenotype","EBV","TGV","TBV","TDD","Residual")
        updateSelectInput(session, inputId = 'ycol', label = 'Y Variable',choices = variables, selected = variables[1])
        updateSliderInput(session,inputId = "startgeneration",label = "Minimum Plot Generation:",min = 0,
                          max = values$df_data[nrow(values$df_data),1],value = 0)
        updateSliderInput(session,inputId = "endgeneration",label = "Maximum Plot Generation:",min = 0,
                          max = values$df_data[nrow(values$df_data),1],value = values$df_data[nrow(values$df_data),1])
      }
      if(datafile == "Inbreeding")
      {
        fileinputed$filename = "Inbreeding";
        values$df_data <- read.table(input$file1$datapath, header=TRUE, sep = " ",colClasses = c("numeric",rep("character",14)))
        temp <- values$df_data ;
        columns <- c(2:8,10:15)
        for(i in 1:length(columns)){temp[ ,columns[i]] <- as.numeric(matrix(unlist(strsplit(temp[,columns[i]],'[()]')),ncol=2,byrow=TRUE)[,1])}
        temp <- apply(temp,2,as.numeric)
        values$df_data <- data.frame(temp,stringsAsFactors = FALSE);
        ## Update columns that can be plotted ##
        variables <- names(values$df_data)
        variables <- variables[2:length(variables)]
        COLUMNS$df_data <- c("Generation","Pedigree_f","Genomic_f","H1_f","H2_f","ROH_f","Homozygosity","ProportionROH","ExpectedHet","Fitness",
                      "Numb_homozlethal","Numb_hetezlethal","Numb_homozysublethal","Numb_hetezsublethal","LethalEquivalents")
        updateSelectInput(session, inputId = 'ycol', label = 'Y Variable',choices = variables, selected = variables[1])
        updateSliderInput(session,inputId = "startgeneration",label = "Minimum Plot Generation:",min = 0,
                          max = values$df_data[nrow(values$df_data),1],value = 0)
        updateSliderInput(session,inputId = "endgeneration",label = "Maximum Plot Generation:",min = 0,
                          max = values$df_data[nrow(values$df_data),1],value = values$df_data[nrow(values$df_data),1])
      }
      if(datafile == "QTL")
      {
        fileinputed$filename = "QTL";
        values$df_data <- read.table(input$file1$datapath, header=TRUE, sep = " ",colClasses = c(rep("numeric",13)))
        ## Update columns that can be plotted ##
        variables <- names(values$df_data)
        variables <- variables[2:length(variables)]
        COLUMNS$df_data <- c("Generation","QTL_Founder_Start","QTL_Founder_Lost","Mutation_QTL_Total","Mutation_QTL_Lost","Additive Variance","Dominance Variance",
                     "FTL_Founder_Start","FTL_Founder_Lost","Mutation_FTL_Total","Mutation_FTL_Lost","Avg_Haplotypes_Within_Window","ProgenyDiedFitness")
        updateSelectInput(session, inputId = 'ycol', label = 'Y Variable',choices = variables, selected = variables[1])
        updateSliderInput(session,inputId = "startgeneration",label = "Minimum Plot Generation:",min = 0,
                          max = values$df_data[nrow(values$df_data),1],value = 0)
        updateSliderInput(session,inputId = "endgeneration",label = "Maximum Plot Generation:",min = 0,
                          max = values$df_data[nrow(values$df_data),1],value = values$df_data[nrow(values$df_data),1])
      }
      if(datafile == "Decay")
      {
        fileinputed$filename = "Decay";
        df <- read.table(file=input$file1$datapath,sep=" ",header=FALSE)
        Distance <- df[1, ]; df <- df[-c(1), ]; # grab distance 
        # Plot for generations 0 5 10 15 20 
        generations <- unique(c(seq(1,nrow(df),5)),nrow(df));
        #generations <- c(1:nrow(df))
        for(i in 1:length(generations)) 
        { 
          tempdf <- data.frame(cbind(t(Distance),t(df[generations[i], ])),stringsAsFactors = FALSE) 
          tempdf$gen <- generations[i] -1; 
          names(tempdf) <- c("Distance","r2","Generation") 
          if(i == 1){plotdf <- tempdf;rm(tempdf)} 
          if(i > 1){plotdf <- rbind(plotdf,tempdf); rm(tempdf)} 
        } 
        plotdf$Generation <- as.factor(as.numeric(plotdf$Generation)) 
        values$df_data <- data.frame(plotdf,stringsAsFactors = FALSE);
        updateSelectInput(session, inputId = 'ycol', label = 'Y Variable',choices = "", selected = "")
        updateSliderInput(session,inputId = "startgeneration",label = "Minimum Plot Generation:",min = 0,
                          max = (nrow(df)-1),value = 0)
        updateSliderInput(session,inputId = "endgeneration",label = "Maximum Plot Generation:",min = 0,
                          max = (nrow(df)-1),value = (nrow(df)-1))
        
      }
    }
  })
  
  output$caption <- renderText({paste(unlist(strsplit(as.character(paste(input$file1[[1]])),"_")),sep=" ",collapse=" ")})
  
  observeEvent(fileinputed$filename, 
  {
    output$MyPlot <- renderPlot(
    {
      if(!is.null(values$df_data))
      {
        if(fileinputed$filename == "Performance" | fileinputed$filename == "Inbreeding" | fileinputed$filename == "QTL")
        {
          x <- values$df_data
          X <- which(names(x) == input$ycol)
          if(length(X) > 0)
          {
            x <- x[ ,c(1,X)]
            x <- data.frame(x)
            # Update Generations to Plot #
            minvalue <- input$startgeneration;
            maxvalue <- input$endgeneration;
            x <- x[which(x[,1] >= minvalue & x[,1] <= maxvalue), ]
            names(x) <- c("X","Y");
           ggplot(x, aes(x=X,y=Y)) + geom_line(size = 1.5) + 
             ggtitle("\n\nTrend Across Time") + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
             ylab (COLUMNS$df_data[X]) + xlab ("Generation") +
             theme(plot.title = element_text(size = 20,hjust = 0.5),axis.text=element_text(size=16),axis.title=element_text(size=14))
          }
        } else if (fileinputed$filename == "Decay") {
          x <- data.frame(values$df_data,stringsAsFactors = FALSE);
          # Update Generations to Plot #
           minvalue <- input$startgeneration;
           maxvalue <- input$endgeneration;
           x$Generation <- as.numeric(paste(x$Generation))
           x <- x[which(x$Generation >= minvalue & x$Generation <= maxvalue), ]
           x$Generation <- as.factor(as.numeric(x$Generation))
           ggplot(x, aes(x=Distance, y = r2, colour = Generation)) + geom_point(size = 1) +
             labs(title = "\n\nLD Decay", x = "Distance (kb)", y = "Correlation between two SNP") + theme_bw() +
             theme(plot.title = element_text(size = 20,hjust = 0.5),axis.title = element_text(size = 16),legend.position="bottom",axis.text=element_text(size=14))
        } else {}
      }
    })
  })
  
  # Show the first "n" observations ----
  #output$view <- renderTable({
  # head(values$df_data, n = input$obs)
  #})
}

shinyApp(ui = ui, server = server)
