library(rsconnect)
library(base)
library(utils)
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


load('scShinyapp081720.RData')


ui <- fluidPage(
    tags$h1("Single-cell spinal cord gene database"),
    tags$h2('For further details, see:'),
    tags$h3('Niehaus et al (2020) Border-associated macrophages resolve pain hypersensitivity after tissue injury, submitted.'),
    setBackgroundColor("aliceblue"),
    br(),
    sidebarLayout(position="left",sidebarPanel(
    textInput(
      inputId = "search", label = "Gene of interest",
      value='Gad2',
      #btnSearch = icon("search"),
      #btnReset = icon("remove"),
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
      else if (toupper(input$search) %in% rownames(df1)){
      propbar_lattice(gene=input$search,sham=shamprop(),sni=sniprop())
        write(input$search,'inputgenes.txt',append=T)}
      else {
        propbar_lattice(gene=input$search,sham=shamprop(),sni=sniprop())
      }
    })}


shinyApp(ui = ui, server = server)
