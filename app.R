# cuePCR: Intuitive analysis of quantitative real-time RT-PCR data.

# Copyright (C) 2020 Dario Angelo Ricciardi
# dario.a.ricciardi@gmail.com

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>


############################################## cuePCR ##############################################################

library(shiny)
library(shinyjs)
library(ggplot2)
library(readxl)
library(car)
library(plyr)


# Define UI############################################################################################################

ui <- navbarPage(title = div(align = "center",
                             img(src="AppLogo_white.png",
                                 style = "margin-top:-5px; width:100px;")),
                 id = "navbar",
                 position = "fixed-top",
                 theme = "bootstrap.css",
                 collapsible = TRUE,
                 inverse = FALSE,
                 windowTitle = "cuePCR | Intuitive qPCR analysis",
                 tabPanel("Start",
                          fluidRow(
                            column(3, 
                                   div(align = "center",
                                       #h1("Shiny qPCR"),
                                       div(img(src="AppLogo_color.png", class = "img_logo"), style = "margin-top:30px; margin-bottom:50px;"),
                                       h4(style = "padding-bottom: 0px;", "Upload data and get started!"),
                                       fileInput("file1","Select *.xls or *.csv result file(s)", multiple = TRUE, tags$style("width:90%;")),
                                       actionButton("readfiles", "Read selected files!",
                                                    style = "width:70%;", class = "btn_std")
                                   ),#End div
                                   div(align = "center", style = "position: fixed; bottom: 50px; width: 25%;",
                                       p("
                                       Did this tool make your work easier? Consider donating to help me cover the server fees and improve cuePCR.
                                         ", style = "text-align:center; font-size:16px;padding:20px;"),
                                   HTML('
                                             <form action="https://www.paypal.com/cgi-bin/webscr" method="post" target="_top">
                                            <input type="hidden" name="cmd" value="_s-xclick" />
                                            <input type="hidden" name="hosted_button_id" value="HGW3LCBLLJG6W" />
                                            <input type="image" src="DonationButtonBlue.png" border="0" name="submit" title="PayPal - The safer, easier way to pay online!" alt="Donate with PayPal button" />
                                            <img alt="" border="0" src="https://www.paypal.com/en_DE/i/scr/pixel.gif" width="1" height="1" />
                                            </form>
                                             ')),
                                   tags$head(tags$meta(name="description", content="A free and easy to use tool for the analysis of qPCR data."),
                                             tags$meta(name="keywords", content="qPCR,cuepcr,gui,analysis,darioar,shiny,dario,ricciardi,plotting,deltact,ct,foldchange,relative,expression"),
                                             tags$meta(name="author", content="Dario A. Ricciardi"),
                                             tags$meta(name="google-site-verification", content="1Lbv2yxufefT0c_kMUHQJOShU52eIKMZ9vbkYFqnxz8"))
                            ),#End column1
                            column(6,
                                   tags$div(align = "center",
                                     h3("From raw data to publication-ready plots in 5 easy steps!", style = "text-align:center; padding:30px; font-weight:bold;")
                                   ),
                                   div(align = "center", img(src = "Instructions4.png", style = "width:100%;"), style = "padding:50px;"),
                                   p(
                                   tags$b("cuePCR"), 
                                   "combines intuitive analysis of real-time quantitative RT-PCR data with a highly customizable
                                   plotting feature. It is free to use and doesn't require any coding knowledge.
                                   Data from the StepOne Plus (Applied Biosystems) cycler can be uploaded as is. Data from other sources must
                                   be in the format of the
                                     ", downloadLink("example_data", "example data set.", style = "font-size:22px;"), 
                                     "To get started, simply upload your data and follow the instructions presented to you in each step."
                                     , style = "text-align:justify; font-size:22px;padding:50px;"),
                                   p(
                                     "Created by Dario Angelo Ricciardi (dario.a.ricciardi@gmail.com)"
                                     , style = "text-align:left; font-size:18px;padding-left:50px;position: fixed; bottom: 30px;") 
                            ),#End column
                            column(width = 1)
                             )#End FluidRow
                 ),#End TabPanel 1
                 tabPanel("Data Management",
                          fluidRow(
                            column(3,  align = "center",
                                   h3("Build Groups", div("?", class = "toolgroup",
                                                          style = "vertical-align:super;font-size:14px;font-weight:bold;",
                                                          # tags$video(tags$source(src = "tutorial.mp4", type = "video/mp4"),
                                                          #      class ="toolgroupvideo")
                                                          HTML('<video width="800" preload muted autoplay loop class = "toolgroupvideo" src = "tutorial.mp4" type = "video/mp4"></video>')
                                   )),
                                   p(style = "text-align:justify;",
                                     "
                                     Select all samples from one condition in the left drop-down menu and the corresponding controls
                                     in the right drop-down menu. Click the Add Group button to confirm your selection. Repeat this
                                     process for all your conditions until all samples are assigned to groups. All samples in one group will be considered as biological replicates.
                                     Check the targets you want to use as reference in the list below. At least one group and one reference gene are required
                                     for the analysis.
                                     "),
                                   div(align = "center", style = "margin-top:50px;",
                                       h3("Select Reference Genes")
                                   ),#End div
                                   uiOutput("refUI"),
                                   useShinyjs(),
                                   div(align = "center",
                                       h4("OPTIONAL: Supply efficiencies for targets", div("?", class = "tooleff",
                                                                                           style = "vertical-align:super;font-size:14px;font-weight:bold;",
                                                                                           span("To supply your own efficiency values, upload a file (xls, csv, txt) with two columns:
                                                                                       The first column should contain a list of your targets as they appear in your data and the second column
                                                                                       should contain the corresponding efficiency values. You can keep a local file of all your targets. The program will
                                                                                       only use the ones that are required for the current experiment.",
                                                                                                class ="toolefftext"))),
                                       
                                       fileInput("efficiency", label = NULL, placeholder = "Default = 2.00"),
                                       actionButton("AnalButt", "Analyse!", style = "width:60%; margin-bottom:5px", class = "btn_std"),
                                       checkboxInput("paired", "Paired data", value = FALSE),
                                       div(style = "margin:60px;")
                                   )
                            ),
                            column(9, align = "center",
                                   fluidRow(
                                     column(5, align = "center",
                                            uiOutput("sampleUI1")
                                     ),
                                     column(2, align = "center",
                                            div(style= "padding-top:60px;"),
                                            uiOutput("addgroupUI"),
                                            uiOutput("removegroupUI")
                                     ),
                                     column(5, align = "center",
                                            uiOutput("sampleUI2")
                                     )
                                   ),#End FluidRow
                                   fluidRow(
                                     column(12, align = "justify",
                                            div(style = "padding-top:50px", uiOutput("RepGroupsUI"))
                                     )
                                   )#End FluidRow
                            )#End column
                          )#End Fluid Row
                 ),#End TabPanel 2
                 tabPanel("Plotting",
                          fluidRow(
                            column(3, align = "center",
                                   h3("Customize Plot"),
                                   selectInput("by", "Choose which parameters to plot by",
                                               choices = c("Target~Group", "Group~Target")),
                                   selectInput("plotvalue", "Choose which value to plot",
                                               choices = c("Log2 FoldChange" = "log2fc", "Foldchange" = "foldchange", "deltaCq" = "dct")),
                                   # selectInput("geom", "Choose graph type",
                                   #             choices = c("Bar", "Boxplot")),
                                   selectInput("palette", "Choose color palette",
                                               choices = list(
                                                 "Sequential" = c("YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu",
                                                                  "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd", "Oranges",
                                                                  "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues"),
                                                 "Qualitative" = c("Set1", "Set2", "Set3", "Pastel1", "Pastel2", "Paired",
                                                                   "Dark2", "Accent"),
                                                 "Diverging" = c("Spectral", "RdYlGn", "RdYlBu", "RdGy", "RdBu", "PuOr",
                                                                 "PRGn", "PiYG", "BrBG")
                                               ),
                                               selected = "Blues"
                                   ),
                                   div(align = "left", style = "padding-left:20%;",
                                       checkboxInput("pal_dir", "Reverse order of colors"),
                                       checkboxInput("contour", "Draw contours", value = TRUE),
                                       checkboxInput("vline", "Show zero-line", value = TRUE)
                                   ),
                                   sliderInput("width", "Bar width", 0.1, 1, 0.8, step = 0.1),
                                   sliderInput("aspectratio", "Plot aspect ratio", 0.1, 2, step = 0.1, value = 0.5),
                                   downloadButton("dload_plot", "Download Plot", style = "width:70%; margin-top:30px;", class = "dload"),
                                   selectInput("dlformat", "Select download format for plot", choices = c("pdf", "png", "tiff", "jpg"))
                            ),#End column 1
                            column(6,
                                   div(align = "center", style = "height:600px;",
                                       uiOutput("plotUI")
                                   ),#End div
                                   fluidRow(
                                     column(width = 4,
                                            selectInput("legend", "Legend position",
                                                        choices = c("none","left","right","top","bottom"),
                                                        selected = "right"),
                                            textInput("legendtitle", "Legend title", value = "", placeholder = "default"),
                                            checkboxInput("legendbox", "Draw legend box", value = FALSE)
                                     ),#End column 2.1
                                     column(width = 4,
                                            textInput("ylab","Y-axis label", value = "", placeholder = "default"),
                                            textInput("xlab", "X-axis label", value = "", placeholder = "default"),
                                            checkboxInput("sigalpha", "Show significance asterisks", value = TRUE)
                                     ),#End column 2.2
                                     column(width = 4,
                                            sliderInput("textsize", "Text size", 8, 32, value = 14, step = 2),
                                            selectInput("errorbarcolor", "Error bar color",
                                                        choices = c("Black", "Grey20", "Grey50", "Grey80", "White", "Red", "Blue", "Green", "Yellow")),
                                            checkboxInput("errorbaralpha", "Error bars", value = TRUE),
                                            div(style = "margin:40px;")
                                            
                                     )#End column 2.3
                                   )#End fluidRow 2
                            ),#End column 2
                            column(3, align = "center",
                                   tags$b("Change group names"),
                                   uiOutput("GroupNameUI")
                            )#End column 3
                          )#End FluidRow 1
                 ),#End TabPanel 3
                 tabPanel("Raw Data",
                          fluidRow(
                            #verbatimTextOutput("sentinel"),
                            uiOutput("dataUI")
                          )#End FluidRow
                 )#End TabPanel 4
)#End NavBarPage



#DEFINE SERVER LOGIC###############################################################################


#Functions#########################################################################################
data.cleanup <- function(x){
  customdata <- FALSE
  if(grepl("\\.xls.*$", x) == TRUE){
    try(fulldata <- read_excel(x, col_names = FALSE))
    try(fulldata <- read_excel(x, col_names = FALSE, sheet = "Results"))
    if(agrepl("Block Type", fulldata[1,1], ignore.case = TRUE) == FALSE){
      customdata <- TRUE
    }else{customdata <- FALSE}
  }else if(grepl("\\.csv$", x) == TRUE){
    fulldata <- read.csv2(x, header = FALSE)
    if(agrepl("sample", fulldata[1,1], ignore.case = TRUE) == TRUE){
      customdata <- TRUE
    }else{customdata <- FALSE
    fulldata[,10] <- gsub(",",".", fulldata[,10])
    fulldata[,11] <- gsub(",",".", fulldata[,11])
    fulldata[,12] <- gsub(",",".", fulldata[,12])
    }
  }
  if(customdata == TRUE){
    
    colnames(fulldata) <- c("sample", "target", "ct_mean")
    if(grepl("[a-z]", fulldata$ct_mean[1]) == TRUE){
    fulldata <- fulldata[-1,]
    fulldata$ct_mean <- gsub(",",".", fulldata$ct_mean)
    
    }else{}
    
    fulldata$ct_mean <- as.numeric(fulldata$ct_mean)
    fulldata$ct_sd <- 0
    fulldata$group <- ""
    fulldata$control <- FALSE
    fulldata$reference <- FALSE
    fulldata$ct_ref <- 0
    fulldata$efficiency <- 2
    return(fulldata)
  }else if(customdata == FALSE){
    fulldata <- data.frame(lapply(fulldata, as.character), stringsAsFactors=FALSE)
    fulldata <- fulldata[9:104, c(1:4, 10:12)]
    colnames(fulldata) <- c("well", "sample", "target", "task", "ct", "ct_mean", "ct_sd")
    fulldata$ct <- as.numeric(fulldata$ct)
    fulldata$ct_mean <- as.numeric(fulldata$ct_mean)
    fulldata$ct_sd <- as.numeric(fulldata$ct_sd)
    fulldata <- subset(fulldata, task == "UNKNOWN")
    fulldata$group <- ""
    fulldata$control <- FALSE
    fulldata$reference <- FALSE
    fulldata$ct_ref <- 0
    fulldata$efficiency <- 2
    fulldata <- fulldata[!duplicated(fulldata[,c("sample", "target", "ct_mean", "ct_sd")]),]
    return(fulldata)
  }
}


datasum <- function(data = NULL, measurevar, groupvars){
  
  sumdata <- ddply(data, groupvars, .drop = TRUE, .fun = function(x, col){
    c(n = length(x[[col]]),
      mean = mean(x[[col]]),
      sd = sd(x[[col]]))
  },
  measurevar
  )
  #Rename mean column
  sumdata <- rename(sumdata, c("mean" = measurevar))
  #Calculate standard error of mean
  sumdata$se <- sumdata$sd / sqrt(sumdata$n)
  
  return(sumdata)
  
}

server <- function(input, output, session) {

#Define reactive variables##########################################################################  
  
  exdata <- reactiveVal()
  inputdata <- reactiveVal()
  inputdata_groups <- reactiveVal()
  stopper <- reactiveVal(FALSE)
  
  gNo <- reactiveVal(value = 1)
  RepGroups <- reactiveVal(value = list())
  resetInput <- reactiveVal(0)
  ref_test <- reactiveVal()
  efficiency <- reactiveVal()
  result_groups <- reactiveVal()
  stat <- reactiveVal()
  statsum <- reactiveVal()
  
  plot_x <- reactiveVal()
  plot_y <- reactiveVal()
  plot_by <- reactiveVal()
  palette_dir <- reactiveVal()
  vlinealpha <- reactiveVal()
  contouralpha <- reactiveVal()
  legendbox <- reactiveVal()
  legendtitle <- reactiveVal()
  x_title <- reactiveVal()
  y_title <- reactiveVal()
  y_title1 <- reactiveVal()
  errorbaralpha <- reactiveVal()
  sigalpha <- reactiveVal()
  gNames <- reactiveVal()
  gNameLength <- reactiveVal()
  
  sentinel <- reactiveVal()
  
  
  
  
  
#Data Input#########################################################################################
  
  #Read Files
  
  #Read example data and store in variable for download by user
  exdata(read_excel("./Example_Data.xlsx"))
  
  #Check for correct file format and give error
  observeEvent(input$file1, {
    if(all(grepl("\\.csv$|\\.xls$|\\.xlsx$", input$file1$datapath)) == TRUE){
    }else{
      showModal(modalDialog(
        div(align = "center",
            h2("Wrong file format!"),
            hr(),
            h3("Only .xls, .xlsx or .csv files are supported")),
        footer = tagList(modalButton("Understood")),
        easyClose = TRUE
      ))
    }
  })
  
  #Enable Read Files button only when correct file type is uploaded
  observe({
    if(length(input$file1$datapath) >= 1){
      if(all(grepl("\\.csv$|\\.xls$|\\.xlsx$", input$file1$datapath)) == TRUE){
        enable("readfiles")
      }else{
        disable("readfiles")
      }}
    else{
      disable("readfiles")
    }
  })
  
  #Read file and create inputdata dataframe
  observeEvent(input$readfiles, {
    removeModal()
    indata <- do.call(rbind, lapply(input$file1$datapath, function(i){
      temp <- data.cleanup(i)
      #If there is no plate column, add one with the corresponding plate number
      if("plate" %in% colnames(temp) == FALSE){
        temp$plate <- as.character(which(input$file1$datapath == i, arr.ind = TRUE))
      }else{}
      temp
    }))
    
    inputdata(indata)
  })
  
  #Data Management####################################################################################
  
  #Detect high standard deviation and give warning
  observeEvent(input$readfiles, {
    
    #Switch to Data Management tab
    updateNavbarPage(session, "navbar", selected = "Data Management")
    
    stopper(FALSE)
    if(stopper() == FALSE){
      if(any(inputdata()$ct_sd > 0.5) == TRUE){
        showModal(modalDialog(
          div(align = "center",
              h2("High standard deviation between technical replicates detected!"),
              hr()),
          h3("Samples:"),
          h3(paste(unique(inputdata()$sample[which(inputdata()$ct_sd > 0.5)]), collapse = ", ")),
          h3("in wells"),
          h3(paste(inputdata()$well[which(inputdata()$ct_sd > 0.5)], collapse = ", ")),
          footer = tagList(actionButton("removeSD", "Remove from analysis"), actionButton("ignoreSD", "Ignore"))
        ))
      }else{stopper(TRUE)}
    }
  })
  
  #Remove SD modal and update stopper variable
  observeEvent(input$ignoreSD, {
    removeModal()
    stopper(TRUE)
  })
  
  observeEvent(stopper(),{
    if(stopper() == TRUE){
      #Check for multiple plates and ask about sample handling
      if(length(unique(inputdata()$plate)) > 1){
        showModal(modalDialog(
          div(align = "center",
              h2(paste0("Multiple files/plates (", length(unique(inputdata()$plate)), ") detected!")),
              hr(),
              h3("Sample handling:")),
          radioButtons("multiplate", "Do you want to add the plate number to each sample name?",
                       choiceNames = c("Yes (this will make data management more complicated, use only if needed)", "No, keep it simple"),
                       choiceValues = c("TRUE", "FALSE"),
                       selected = "FALSE", width = '90%'),
          footer = tagList(modalButton("Apply"))
        ))
      }
    }else{}
  })
  
  #Remove high SD samples from dataset
  observeEvent(input$removeSD, {
    indata <- inputdata()
    indata <- indata[-which(indata$ct_sd > 0.5),]
    inputdata(indata)
    removeModal()
    stopper(TRUE)
  })
  
  #Add plate numbers for multiplate data
  observeEvent(input$multiplate, {
    if(is.null(input$multiplate) == FALSE){
      if(input$multiplate == TRUE){
        indata <- inputdata()
        indata$sample <- paste0(indata$sample, "@plate", indata$plate)
        inputdata(indata)
      }else{}
    }else{}
  })
  
  #Render sample UI: Select list input 1
  output$sampleUI1 <- renderUI({
    if(length(input$file1$datapath) >= 1){
      tagList(
        h3("Condition Samples"),
        selectInput("condition", "Select all samples from one condition", choices = unique(inputdata()$sample), multiple = TRUE, tags$style("width:80%;"))
      )
    }else{}
  })
  
  #Render sample UI: Add group button
  output$addgroupUI <- renderUI({
    if(length(input$file1$datapath) >= 1){
      actionButton("addgroup", h3(tags$b("Add Group")), style = "width:180px;", class = "btn_add")
    }else{
      "Upload a file to see samples"
    }
  })
  
  #Render sample UI: Select list input 2
  output$sampleUI2 <- renderUI({
    if(length(input$file1$datapath) >= 1){
      tagList(
        h3("Control Samples"),
        selectInput("control", "Select corresponding control sample(s)", choices = unique(inputdata()$sample), multiple = TRUE, tags$style("width:80%;"))
      )
    }else{}
  })
  
  #Render Replicate Group UI: Group tiles
  output$RepGroupsUI <- renderUI({
    
    if(length(input$file1$datapath) >= 1){
      div(align = "center",
          lapply(RepGroups(), function(i){
            actionButton("placeholder", label = 
                           div(
                             p(paste0("Group ",unique(i$Group)), style = "margin-top:0px; margin-bottom:-5px; font-size:18px;"),
                             hr(style = "margin-bottom: 10px; margin-top:-5px;"),
                             p(paste(unique(i$Condition), collapse = ", "), style = "margin:0px;"),
                             p("vs.", style = "margin:0px"),
                             p(paste(unique(i$Control), collapse = ", "), style = "margin:0px;"), style = "padding:0px;"
                           ), class = "btn_group")
          })
      )
    }else{
      ""
    }
    
  })
  
  #Add Group of selected samples to RepGroups
  observeEvent(input$addgroup,{
    if(is.null(input$control) == TRUE){
      
    }else if(is.null(input$control) == FALSE){
      if(length(input$control) != 0 & length(input$condition) != 0){
        if(any(input$condition %in% input$control) == TRUE){
          showModal(modalDialog(
            div(align = "center",
                h2("Can not compare a sample to itself!"),
                hr(),
                h3("Please choose different samples for condition and control")),
            footer = tagList(modalButton("Understood")),
            easyClose = TRUE
          ))
          resetInput(resetInput()+1)
        }else{
          temp <- RepGroups()
          temp[[gNo()]] <- list(Group = gNo(), Control = input$control, Condition = input$condition)
          RepGroups(temp)
          gNo(gNo()+1)
        }
      }
    }
  })
  
  #Reset select inputs to blank
  observeEvent({RepGroups()
    resetInput()},{
      updateSelectInput(session, "control", selected = character(0))
      updateSelectInput(session, "condition", selected = character(0))
    })
  
  
  #Reset RepGroups and group number counter (gNo) upon file selection
  observeEvent(input$file1,{
    RepGroups(list())
    gNo(1)
  })
  
  #Remove last group from RepGroups and decrease gNo
  observeEvent(input$removegroup,{
    gNo(gNo()-1)
    temp <- RepGroups()
    temp[[gNo()]] <- NULL
    RepGroups(temp)
  })
  
  #Render Remove group button
  output$removegroupUI <- renderUI({
    if(gNo() > 1){
      actionButton("removegroup", "Remove Group", style = "width:70%;", class = "btn_rm")
    }else{}
  })
  
  #Render RepGroups as data table by binding individual entries of RepGroups list
  output$RepGroups <- renderTable({
    do.call(rbind, RepGroups())
  })
  
  #Disable Analyse button if no Groups have been added yet and/or no reference gene has been picked yet
  observe({
    if(gNo() > 1 & any(ref_test()) == TRUE){
      enable("AnalButt")
    }else{disable("AnalButt")}
  })
  
  
  #Render Reference gene selection
  output$refUI <- renderUI({
    if(length(input$file1$datapath) == 0){
      div(style = "margin:10px 0px 50px 0px;",
          "Upload a file to see reference genes")
    }else{
      div(align = "left", style = "margin:20px 0px 40px 25%;",
          lapply(unique(inputdata()$target), function(i) {
            checkboxInput(inputId = paste0("target_", i),
                          label = paste0(i)
            )
          })
      )
    }
  })
  
  
  #Monitor Reference gene selection status and enable Analyse button only if at least one gene is selected
  observe({
    ref_test1 <- ref_test()
    for(i in unique(inputdata()$target)){
      if(isTRUE(input[[paste0("target_", i)]]) ==  TRUE){
        ref_test1[[i]] <- TRUE
      }else{
        ref_test1[[i]] <- FALSE
      }
    }
    ref_test(ref_test1)
  })
  
  #Read efficiency file
  observeEvent(input$efficiency,{
    if(grepl("\\.xls.*$", input$efficiency$datapath) == TRUE){
      eff <- read_excel(input$efficiency$datapath, col_names = FALSE)
      colnames(eff) <- c("target", "efficiency")
      efficiency(eff)
    }else if(grepl("\\.csv$", input$efficiency$datapath) == TRUE){
      eff <- read.csv2(input$efficiency$datapath, header = FALSE)
      eff[,2] <- gsub(",",".", eff[,2])
      eff[,2] <- as.numeric(eff[,2])
      colnames(eff) <- c("target", "efficiency")
      efficiency(eff)
    }else if(grepl("\\.txt$", input$efficiency$datapath) == TRUE){
      eff <- read.table(input$efficiency$datapath)
      eff[,2] <- gsub(",",".", eff[,2])
      eff[,2] <- as.numeric(eff[,2])
      colnames(eff) <- c("target", "efficiency")
      efficiency(eff)
    }else{
      showModal(
        modalDialog({
          div(align = "center",
              h3("Invalid file type. Allowed formats: xls, xlsx, csv, txt"))
        })
      )
    }
    
  })
  
  #Calculation after pressing Analyse Button
  observeEvent(input$AnalButt, {
    #Prepare data
    indata <- inputdata()
    indata$group <- ""
    indata$Group <- ""
    indata$control <- FALSE
    
    #Insert efficiencies
    eff_all <- efficiency()
    eff <- eff_all[eff_all$target %in% unique(indata$target),]
    for (i in unique(eff$target)) {
      indata$efficiency[indata$target %in% i] <- eff$efficiency[eff$target %in% i]
    }
    
    #Assign samples to previously determined groups
    for (i in 1:length(RepGroups())) {
      indata$group[indata$sample %in% RepGroups()[[i]]$Condition] <- paste(indata$group[indata$sample %in% RepGroups()[[i]]$Condition] ,RepGroups()[[i]]$Group, sep = ",")
      indata$group[indata$sample %in% RepGroups()[[i]]$Control] <- paste(indata$group[indata$sample %in% RepGroups()[[i]]$Control] ,RepGroups()[[i]]$Group, sep = ",")
    }
    
    #Remove leading comma. Artifact of multiple group assignment
    indata$group <- gsub("^,","", indata$group)
    
    #Give warning about unassigned samples
    if(any(indata$group == "") == TRUE){
      showModal(modalDialog(
        h2("These samples are not assigned to a group and will be excluded from analysis!"),
        hr(),
        h4(paste(unique(indata$sample[which(indata$group == "")]), collapse = ", "))
      ))}else{}
    
    #Mark references
    for(i in unique(indata$target)){
      if(isTRUE(input[[paste0("target_", i)]]) ==  TRUE){
        indata$reference[indata$target %in% i] <- TRUE
      }else{
        indata$reference[indata$target %in% i] <- FALSE
      }
    }
    
    #Calculate mean reference ct per sample and write to new column
    indata_ref <- subset(indata, reference == TRUE)
    for (i in unique(indata$sample)) {
      indata$ct_ref[indata$sample %in% i] <- mean(indata_ref$ct_mean[indata_ref$sample %in% i], na.rm = TRUE)
    }
    
    #Calculate dct (ct_target - ct_reference)
    indata$dct <- indata$ct_mean - indata$ct_ref
    
    
    #Quick Maffs!
    temp_indata <- data.frame()
    stat <- data.frame("group" = character(), "Target" = character(), "test" = character(), "pval" = double())
    
    for (i in 1:length(RepGroups())) {
      temp <- subset(indata, grepl(paste0("\\b",i,"\\b"), indata$group)) #Set word boundaries (\\b) for grep pattern to get the exact group number
      temp$dctRQ <- 0
      temp$RQ <- 0
      temp$foldchange <- 0
      temp$log2fc <- 0
      #Mark controls
      temp$control[temp$sample %in% RepGroups()[[i]]$Control] <- TRUE
      
      
      #statistics between Control and Condition deltaCt
      tempnoref <- subset(temp, reference == FALSE)
      
      for (t in unique(tempnoref$target)) {
        try({
          temptarget <- subset(tempnoref, target == t)
          tempctrl <- subset(temptarget, control == TRUE)
          tempcond <- subset(temptarget, control == FALSE)
          sentinel(tempctrl)
          if(length(na.omit(tempctrl$dct)) >= 3 & length(na.omit(tempcond$dct >= 3))){
            #test for normality
            #if(shapiro.test(tempctrl$dct)$p.value > 0.05 & shapiro.test(tempcond$dct)$p.value > 0.05){
            #normal distribution -> t-Test
            ttest <- t.test(tempctrl$dct,
                            tempcond$dct,
                            var.equal = if(leveneTest(data = tempnoref, dct ~ control)$`Pr(>F)`[1] < 0.05){FALSE}else{TRUE},
                            paired = input$paired)
            
            tempframe <- data.frame("group" = i, "Target" = t, "test" = ttest$method, "pval" = ttest$p.value)
            stat <- rbind(stat, tempframe)
          }else{
            tempframe <- data.frame("group" = i, "Target" = t, "test" = "Not enough replicates (>= 3 required)", "pval" = NA)
            stat <- rbind(stat, tempframe)
          }       
          # }else{
          #   #no normal distribution -> Wilcoxon Mann Whitney U Test
          #   whitneyhouston <- wilcox.test(tempctrl$dct, tempcond$dct, paired = input$paired)
          #   tempframe <- data.frame("group" = i, "Target" = t, "test" = whitneyhouston$method, "pval" = whitneyhouston$p.value)
          #   stat <- rbind(stat, tempframe)
          # }
        })
      }
      
      #Split into ctrl and condition group
      temp_ctrl <- subset(temp, control == TRUE)
      temp_cond <- subset(temp, control == FALSE)
      
      #Calculate dctRQ for each Target
      for (t in unique(temp_cond$target)) {
        temp_cond$dctRQ[temp_cond$target %in% t] <- mean(temp_ctrl$ct_mean[temp_ctrl$target %in% t], na.rm = TRUE) - temp_cond$ct_mean[temp_cond$target %in% t]
      }
      
      #Split each group into GOI and REF
      temp_GOI_cond <- subset(temp_cond, reference == FALSE)
      temp_REF_cond <- subset(temp_cond, reference == TRUE)
      
      #Calculate RQ values
      temp_GOI_cond$RQ <- temp_GOI_cond$efficiency^temp_GOI_cond$dctRQ
      temp_REF_cond$RQ <- temp_REF_cond$efficiency^temp_REF_cond$dctRQ
      
      #Calculate ratio (foldchange)
      temp_GOI_cond$foldchange <- temp_GOI_cond$RQ / mean(temp_REF_cond$RQ, na.rm = TRUE)
      
      #Calculate log2 foldchange
      temp_GOI_cond$log2fc <- log2(temp_GOI_cond$foldchange)
      
      temp_GOI_cond$group <- i
      temp_REF_cond$group <- i
      temp_ctrl$group <- i
      
      temp_indata <- rbind(temp_indata, temp_GOI_cond, temp_REF_cond, temp_ctrl)
    }
    
    
    #Reassign reactive variables and decouple result dataframes
    stat(stat)
    result_groups(temp_indata)
    inputdata_groups(indata)
    
    #Switch to Plotting tab
    updateNavbarPage(session, "navbar", selected = "Plotting")
  })
  
  
  #Plotting###########################################################################################  
  
  
  #Change plotting formula
  observe({
    if(input$by == "Group~Target"){
      plot_x("Group")
      plot_by("Target")
    } else if(input$by == "Target~Group"){
      plot_x("Target")
      plot_by("Group")
    }
  })  
  
  #Change color palette direction
  observe({
    if(input$pal_dir == FALSE){
      palette_dir(1)
    }else if(input$pal_dir == TRUE){
      palette_dir(-1)
    }
  })
  
  #Add zero line
  observe({
    if(input$vline == TRUE){
      vlinealpha(1)
    }else if(input$vline == FALSE){
      vlinealpha(0)
    }
  })  
  
  #Add contours
  observe({
    if(input$contour == TRUE){
      contouralpha(1)
    }else if(input$contour == FALSE){
      contouralpha(0)
    }
  })  
  
  #Add legend box
  observe({
    if(input$legendbox == TRUE){
      legendbox("solid")
    }else if(input$legendbox == FALSE){
      legendbox("blank")
    }
  })
  
  #Add X-axis title
  observe({
    if(input$xlab == ""){
      x_title(plot_x())
    }else{
      x_title(input$xlab)
    }
  })
  
  #Assign Y-axis titles from plotvalue
  observe({
    if(input$plotvalue == "log2fc"){
      y_title1("Log2 Foldchange")
    }else if(input$plotvalue == "foldchange"){
      y_title1("Foldchange")
    }else if(input$plotvalue == "dct"){
      y_title1(expression(Delta*"Cq"))
    }
  })
  
  #Add Y-axis title
  observe({
    if(input$ylab == ""){
      y_title(y_title1())
    }else{
      y_title(input$ylab)
    }
  })
  
  #Add legend title
  observe({
    if(input$legendtitle == ""){
      legendtitle(plot_by())
    }else{
      legendtitle(input$legendtitle)
    }
  })
  
  #Error Bars
  observe({
    if(input$errorbaralpha == TRUE){
      errorbaralpha(1)
    }else{
      errorbaralpha(0)
    }
  })
  
  #Significance Indicators
  observe({
    if(input$sigalpha == TRUE & input$plotvalue != "dct"){
      sigalpha(1)
    }else{
      sigalpha(0)
    }
  })
  
  #Monitor number of groups for group name inputs
  observe({
    gNameLength(length(unique(result_groups()$group)))
  })
  
  #Render group name inputs dynamically
  output$GroupNameUI <- renderUI({
    if(gNameLength() >= 1){
      isolate({
        #read sample names for each group
        groupnames <- lapply(1:gNameLength(), function(i){
          result_groups()$sample[result_groups()$group == i][1]
        })
        #fill group name inputs with sample names, check for duplicates and add index
        lapply(1:gNameLength(), function(i){
          textInput(inputId = paste0("group_", i),
                    label = h5(paste("Group", i)),
                    value = if(any(duplicated(groupnames)) == TRUE){paste0(groupnames[i], "_", i)
                    }else{groupnames[i]},
                    placeholder = "assign name")
        })
      })
    }else{}
  })
  
  #Read group name inputs dynamically and assign new names to data
  observe({
    gNames <- list()
    result_groups <- result_groups()
    
    gNames <- lapply(unique(result_groups$group), function(i){
      input[[paste0("group_", i)]]
    })
    
    for (i in 1:length(gNames)) {
      result_groups$Group[which(result_groups$group == i)] <- as.character(gNames[i])
    }
    
    result_groups(result_groups)
    gNames(gNames)
  })
  
  
  #Render plot
  output$plot <- renderPlot({
    if(length(result_groups()) > 1){
      
      #Read data for manipulation and plotting
      plotdata <- subset(result_groups(), reference == FALSE & control == FALSE)
      #Rename Target column
      if(colnames(plotdata)[3] == "target"){
      colnames(plotdata)[3] <- "Target"
      }else{colnames(plotdata)[2] <- "Target"}
      #sort group names according to group number instead of alphabetically
      plotdata$Group <- factor(plotdata$Group, levels = unique(plotdata$Group))
      #Summarize plotdata for errorbars
      sumdata <- datasum(subset(plotdata, reference == FALSE & control == FALSE),
                         measurevar = input$plotvalue,
                         groupvars = c("Target", "group", "Group"))
      #Build dataframe for significance indicators
      sumdata <- join(sumdata, stat(), type = "left", by = c("group", "Target"))
      sumdata$sig <- ""
      sumdata$sig[sumdata$pval <= 0.05] <- "*"
      statsum(sumdata)
      
      #Bar plot
      #if(input$geom == "Bar"){
      ggplot(data = sumdata, aes(x = !!as.name(plot_x()),
                                 y = !!as.name(input$plotvalue),
                                 fill = !!as.name(plot_by())))+
        geom_bar(stat = "identity", position = position_dodge(), width = input$width, color = alpha(colour = "black", alpha = contouralpha()))+
        geom_hline(yintercept = 0, color = "black", alpha = vlinealpha())+  
        scale_fill_brewer(palette = input$palette, direction = palette_dir())+
        geom_errorbar(aes(ymin = !!as.name(input$plotvalue)-se, ymax = !!as.name(input$plotvalue)+se),
                      width = input$width/5,
                      position = position_dodge(input$width), color = alpha(input$errorbarcolor, alpha = errorbaralpha()))+
        geom_text(aes(x = !!as.name(plot_x()),
                      y = (abs(!!as.name(input$plotvalue)) + se + 0.3) / abs(!!as.name(input$plotvalue)) * !!as.name(input$plotvalue),
                      group = !!as.name(plot_by()),
                      label = sig),
                  position = position_dodge(width = input$width),
                  color = alpha(colour = "black", alpha = sigalpha()),
                  size = 7)+
        theme_minimal()+
        theme(legend.position = input$legend,
              legend.background = element_rect(fill="white", size=.5, linetype = legendbox()),
              text = element_text(size = input$textsize))+
        guides(fill=guide_legend(title = legendtitle()))+
        ylab(y_title())+
        xlab(x_title())+
        coord_fixed(ratio = input$aspectratio)
      
      #Box plot  
      # }else if(input$geom == "Boxplot"){
      #   ggplot(data = plotdata, aes(x = !!as.name(plot_x()),
      #                             y = !!as.name(input$plotvalue),
      #                             fill = !!as.name(plot_by())))+
      #     geom_boxplot(stat = "identity", width = input$width)+
      #       geom_hline(yintercept = 0, color = "black", alpha = vlinealpha())+  
      #       scale_fill_brewer(palette = input$palette, direction = palette_dir())+
      #      theme_minimal()+
      #       theme(legend.position = input$legend,
      #             legend.background = element_rect(fill="white", size=.5, linetype = legendbox()),
      #             text = element_text(size = input$textsize))+
      #       guides(fill=guide_legend(title = legendtitle()))+
      #       ylab(y_title())+
      #       xlab(x_title())+
      #       coord_fixed(ratio = input$aspectratio)
      # }
    }else{}
  })
  
  #Render Plot UI
  output$plotUI <- renderUI(plotOutput("plot", width = 900, height = 600)
  )
  
  #Download Example Data
  output$example_data <- downloadHandler(
    filename = "Example_Data.csv",
    content = function(file) {
      write.csv2(exdata(), file, row.names = FALSE, sep = "\t", quote = FALSE)}
  )
  
  #Download Result Data
  output$dload_data <- downloadHandler(filename = function(){
    paste0(gsub("\\..*$", "", input$file1$name), "_DATA", ".csv")
  },
  content = function(file) {
    write.csv2(result_groups(), file, row.names = FALSE, sep = "\t", quote = FALSE)}
  )
  
  #Download Statistical Data
  output$dload_stat <- downloadHandler(filename = function(){
    paste0(gsub("\\..*$", "", input$file1$name), "_STAT", ".csv")
  },
  content = function(file) {
    write.csv2(statsum(), file, row.names = FALSE, sep = "\t", quote = FALSE)}
  )
  
  #Download Plot
  output$dload_plot <- downloadHandler(filename = function(){
    paste0(gsub("\\..*$", "", input$file1$name), "_PLOT.", input$dlformat)
  },
  content = function(file) {
    ggsave(file, device = input$dlformat, width = 30, height = 30, units = "cm")}
  )
  
  #Raw Data###########################################################################################
  
  output$dataUI <- renderUI(
    if(length(result_groups()) > 1){
      tagList(
        downloadButton("dload_data", "Download result data", style = "margin:10px"),
        tableOutput("resultgroups"),
        downloadButton("dload_stat", "Download statistical summary", style = "margin:10px"),
        tableOutput("statsum"),
        div(style = "padding:50px")
      )
    }else{}
  )
  
  output$resultgroups <- renderTable({
    result_groups()
  })
  
  output$statsum <- renderTable({
    statsum()
  })
  
  # output$sentinel <- renderPrint({
  #   sentinel()
  # })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
