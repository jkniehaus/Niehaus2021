#online app; first only averages
#install.packages(c('gtools','grid','gtable','ggplot2','shiny','shinyWidgets'), repos='https://cloud.r-project.org')
#setwd("~/Documents/OneDrive/OneDrive - University of North Carolina at Chapel Hill/Paper/Online_Resource")
library(rsconnect)
library(gtools)
library(grid)
library(gtable)
library(ggplot2)
library(shiny)
library(shinyWidgets)
library(data.table)
library(lattice)
library(devtools)
library(trelliscope)
library(ggplotify)
library(profvis)
library(gridBase)
library(gridExtra)


load('scShinyapp.RData')
#load('thisisbullshit2.RData')
#print('latest bullshit2')
#print(dim(shamavg))
#app code
## Not run: 
#if (interactive()) {

ui <- fluidPage(
    tags$h1("Single-cell spinal cord gene database"),
    setBackgroundColor("aliceblue"),
    br(),
    sidebarLayout(position="left",sidebarPanel(
    searchInput(
      inputId = "search", label = "Gene of interest (case sensitive)",
      placeholder = "",
      btnSearch = icon("search"),
      btnReset = icon("remove"),
      width = "450px"
    )
    ),mainPanel(
      fluidRow(splitLayout(cellWidths=c("50%", "50%"),
                           plotOutput('P1',height=900,width=450), plotOutput('P2',height=900,width=450))
    )
    ))
)
server <- function(input, output, session) {
    #df1=read.table('ShamClusterAverages.txt',header=T,sep='\t',row.names=1)
    df1=data.frame(fread('ShamClusterAverages.txt',header=T))
    rownames(df1)=df1$V1
    df1$V1=NULL
    #df2=read.table('SNIClusterAverages.txt',header=T,sep='\t',row.names=1)
    df2=data.frame(fread('SNIClusterAverages.txt',header=T))
    rownames(df2)=df2$V1
    df2$V1=NULL
    shamavg=reactive(df1)
    sniavg=reactive(df2)
    #df3=read.table('ShamProportions.txt',header=T,sep='\t',row.names=1)
    #df4=read.table('SNIProportions.txt',header=T,sep='\t',row.names=1)
    df3=data.frame(fread('ShamProportions.txt',header=T))
    rownames(df3)=df3$V1
    df3$V1=NULL
    df3=df3*100
    df4=data.frame(fread('SNIProportions.txt',header=T))
    rownames(df4)=df4$V1
    df4$V1=NULL
    df4=df4*100
    shamprop=reactive(df3)
    sniprop=reactive(df4)
      

      output$P1 <- renderPlot({
        input$search
        if(input$search==""){}
        else{
        
      avgbar_lattice(gene=input$search,sham=shamavg(),sni=sniavg())
    }})
    output$P2 <- renderPlot({
      input$search
      if(input$search==""){}
      else{
      propbar_lattice(gene=input$search,sham=shamprop(),sni=sniprop())}
    })}


shinyApp(ui = ui, server = server)


## End(Not run)
