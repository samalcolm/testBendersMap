library(shiny)
library(shinyWidgets)
library(dplyr)
library(purrr)
library(igraph)
library(ggplot2)
library(ggraph)
library(ggmap)

nodefilename <- "./nodes.ben"
putfilename <- "./iterations.out" 
convfilename <- "./convergence.out"
prev_time <- 0

usa_map <- geom_polygon(data = map_data('usa'), aes(x = long, y = lat, group = group),
                               fill = "#CFCFCF", color = "#505050", size = 0.15)

maptheme <- theme(panel.grid = element_blank()) +
  theme(axis.text = element_blank()) + 
  theme(axis.ticks = element_blank()) + 
  theme(axis.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_rect(fill = "lightblue")) 

nodes <- read.csv(nodefilename, header = TRUE,
                  quote = "'", sep = "",
                  col.names = c('id', 'name', 'lon', 'lat', 'Category'))

node_map <- geom_point(data=nodes, aes(x = lon, y = lat, fill=Category),                        
           shape = 21, size = 3,
           color = 'black', stroke = 0.5) 

node_labels <- geom_text(data=nodes, aes(x = lon, y = lat, label = name),             
             hjust = 0, nudge_x = 1, #nudge_y = 1, 
             size = 3, color = "black", fontface = "bold")  
  
ui <- fluidPage(
  titlePanel("Benders decomposition progress"),
  fluidRow( 
    column(4, offset=4,
     sliderInput("iterations",h3("Iterations"),value=0,min=0,max=1,step=1,ticks=FALSE, animate=T)
     )
    ),
  fluidRow( 
    wellPanel(
     plotOutput("solutionPlot", width="720px",height="480px")
     )
    ),
  fluidRow( 
    wellPanel(  
     plotOutput("convergencePlot", width="720px",height="300px")
     )
    )
  )

server <- function(input, output, session) {
  # Don't activate slider until finished???
# Note: If final solution occurs with others during an update, it won't be reflected in plot
  updateperiod <- 250 # milliseconds
  solutionData <- reactiveFileReader(updateperiod, session, putfilename, readLines)
  convergenceData <- reactiveFileReader(updateperiod, session, convfilename, readLines)
  
  output$solutionPlot <- renderPlot({
    # Read the solution table as it is updated
    # For larger tables, GAMS could only output current iteration and we would append it here
    text <- solutionData() 
    data <- read.delim(text = text, header = T,
               quote = "'", sep = ",",
               col.names = c('iter', 'from', 'to', 'Quantity'), strip.white=T)
    iterlist <- unique(data$iter)
    
    # while GAMS is still running (no Benders convergence), update slider to mark progress
    if (file.info(putfilename)$mtime != prev_time) {
      updateSliderInput(session, "iterations", value=length(iterlist), max=length(iterlist))
    } 
    
    prev_time <<- file.info(putfilename)$mtime
    
    s2 <- merge(data,nodes[,c("id","lon","lat")],by.x="from",by.y="id")
    names(s2)[5] <- "x"
    names(s2)[6] <- "y"
    s3 <- merge(s2,nodes[,c("id","lon","lat")],by.x="to",by.y="id")
    names(s3)[7] <- "xend"
    names(s3)[8] <- "yend"
    
    p <- ggplot(nodes) + usa_map +
      geom_curve(aes(x = x, y = y, xend = xend, yend = yend, size = Quantity),
                 data = s3[s3$iter==iterlist[input$iterations],], 
                 color="green4",
                 curvature = -0.2, 
                 alpha = 0.75,
                 arrow = arrow(length = unit(0.03, "npc")) )+
      scale_size_continuous( range = c(0.25, 2)) +
      node_map + node_labels +
      maptheme

    return(p)
  }, width=720, height=480)

  output$convergencePlot <- renderPlot({
    text <- convergenceData() 
    data <- read.delim(text = text, header = T,
                       quote = "'", sep = ",",
                       col.names = c('iter', 'UB', 'LB', 'new'), strip.white=T)
    convdata <- reshape2::melt(data,id=c("iter","new"),value.var=c("UB","LB"))
    
    p <- ggplot(convdata) + 
      geom_line(aes(iter,value, group=variable,color=variable)) +
      geom_point(data=convdata[convdata$new==1,],aes(iter,value)) +
      xlab("Iteration") + ylab("Bounds") + theme(legend.title=element_blank())
    
    return(p)
  }, width=720, height=300)
  
}

# Run the application 
shinyApp(ui = ui, server = server)

