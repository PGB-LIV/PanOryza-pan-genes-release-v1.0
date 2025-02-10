library(shiny)
library(ComplexHeatmap)
library(shinyjs)

source("helper.R")

ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  titlePanel("PanOryza Interactive Heatmaps"),
  sidebarLayout(
    sidebarPanel(
      selectInput("heatmap", "Select chromosome:", choices = paste("Chromosome", 1:12))
    ),
    mainPanel(
      actionButton("open_heatmap", "Open Heatmap"),
      h4("Instructions to view genomic position of a pangene:"),
      p("Select a chromosome from the dropdown menu and click 'Open Heatmap' to open it in a new tab."),
      p("Note: Only Pan-genes with occupancy > 2 can be searched and displayed"),
      plotOutput("heatmap_plot", width = "100%", height = "100%")
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$open_heatmap, {
    selected_heatmap <- as.numeric(gsub("Chromosome ", "", input$heatmap))
    ht <- draw(plot_by_chr[[selected_heatmap]])
    
    # Open heatmap in a modal dialog
    showModal(modalDialog(
      title = paste("Heatmap of chromosome", selected_heatmap),
      htShiny(ht),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),
      tags$script(HTML('
        $(document).on("shown.bs.modal", function() {
          // Target the iframe and set its dimensions
          var iframe = $(".shiny-frame");
          iframe.css({
            "width": "calc(100vw - 2px)",
            "height": "calc(100vh - 2px)"  
          });
          $(".modal-dialog").css({
            "width": "100%",
            "height": "100%",
            "margin": "0",
            "padding": "0"
          });
          $(".modal-content").css({
            "height": "100%",
            "border-radius": "0"
          });
          $(".modal-body").css({
            "height": "calc(100% - 55px)",
            "overflow-y": "auto"
          });
        });
      '))
    ))
  })
}

shinyApp(ui, server )
