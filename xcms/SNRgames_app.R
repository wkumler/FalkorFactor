#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)


# Define UI ----
ui <- fluidPage(
   
   # Application title
   titlePanel("Signal-to-noise ratio games"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        sliderInput("noise_mean",
                    "Log10 background intensity",
                    min = 0,
                    max = 6,
                    value = 3,
                    step = 0.1),
        textInput("ylim",
                  label = "Y-axis maximum", 
                  value="Unfixed",
                  placeholder = "Usually = peak height"),
        sliderInput("noise_sd",
                    "Log10 background noise",
                    min = 0,
                    max = 6,
                    value = 2,
                    step = 0.1),
         sliderInput("peak_height",
                     "Log10 scale peak height (above average background)",
                     min = 2,
                     max = 6,
                     value = 3,
                     step = 0.1),
         sliderInput("peak_width",
                     "Peak width?",
                     min = 5,
                     max = 100,
                     value = 50),
        radioButtons("peak_shape",
                     "Peak shape?",
                     choices = list("Gaussian", "Tailed", "Flat")),
         sliderInput("scan_length",
                     "Scan length?",
                     min = 100,
                     max = 1000,
                     value = 500)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        fluidRow(
          column(width = 12,
                 htmlOutput("sideinfo")
          ),
          column(width = 6,
                 plotOutput("distPlot")),
          column(width = 6,
                 plotOutput("xcmsPlot"))
          )
        )
      )
)

# Define server logic ----
server <- function(input, output) {
  output$sideinfo <- renderText({
    HTML(paste("Background intensity:", floor(10^input$noise_mean)),
         paste("Peak height:", floor(10^input$peak_height)),
         paste("Peak width:", input$peak_width),
         paste("Scan length:", input$scan_length),
         paste("Noise SD:", floor(10^input$noise_sd)))
  })

  peak_values <- eventReactive({
    input$peak_width
    input$peak_height
    input$peak_shape
    input$scan_length
    input$noise_sd
    input$noise_mean
  }, {
    if(input$peak_shape=="Gaussian"){
      peak_values <- rnorm(floor(10^input$peak_height)/4*input$peak_width)
    } else if(input$peak_shape=="Tailed"){
      peak_values <- rbeta(floor(10^input$peak_height)/3*input$peak_width, 2, 10)
    } else if(input$peak_shape=="Flat"){
      peak_values <- runif(floor(10^input$peak_height)*input$peak_width)
    }
    
    peak_values <- c(rep(0, (input$scan_length-input$peak_width)/2), 
                     table(cut(peak_values, breaks = input$peak_width)), 
                     rep(0, (input$scan_length-input$peak_width)/2)) %>%
      `+`(rnorm(n = input$scan_length, sd = floor(10^input$noise_sd), mean = floor(10^input$noise_mean)))
    return(peak_values)
    print(peak_values)
  })
  
  output$distPlot <- renderPlot({
    peak_start_index <- input$scan_length/2-input$peak_width/2
    peak_end_index <- input$scan_length/2+input$peak_width/2
    noise_window_indices <- list(c(peak_start_index-input$peak_width,
                                   peak_start_index),
                                 c(peak_end_index,
                                   peak_end_index+input$peak_width))

    left_noise_window_indices <- noise_window_indices[[1]][1]:noise_window_indices[[1]][2]
    left_noise_window_indices <- left_noise_window_indices[left_noise_window_indices>0]
    right_noise_window_indices <- noise_window_indices[[2]][1]:noise_window_indices[[2]][2]
    right_noise_window_indices <- right_noise_window_indices[right_noise_window_indices<input$scan_length]


    sd_noise_estimate <- sd(c(peak_values()[left_noise_window_indices],
                              peak_values()[right_noise_window_indices]),
                            na.rm = T)
    median_noise_estimate <- median(c(peak_values()[left_noise_window_indices],
                                      peak_values()[right_noise_window_indices]),
                                    na.rm = T)

    if(input$ylim=="Unfixed"){
      plot(peak_values(), pch=19)
    } else {
      plot(peak_values(), pch=19, ylim=c(0, as.numeric(input$ylim)))
    }

    abline(h=max(peak_values()), col = "blue")
    mtext(text = "Max\npeakheight", side = 4, at = max(peak_values()),
          col = "blue", line = 1)
    arrows(x0 = noise_window_indices[[1]][1], x1 = noise_window_indices[[1]][2],
           y0 = median_noise_estimate, y1 = median_noise_estimate,
           col = "blue", angle = 90, code = 3, lwd = 2)
    arrows(x0 = noise_window_indices[[2]][1], x1 = noise_window_indices[[2]][2],
           y0 = median_noise_estimate, y1 = median_noise_estimate,
           col = "blue", angle = 90, code = 3, lwd = 2)
    text("Noise window", x = noise_window_indices[[2]][2],
         y = median_noise_estimate+floor(10^input$peak_height)/10, col = "blue")

    legend("topleft", legend = c("Signal-to-noise ratio",
                                 paste("Max/Median:", round(max(peak_values())/median_noise_estimate)),
                                 paste("Max/SD:", round(max(peak_values())/sd_noise_estimate)),
                                 paste("Max-Mean/SD:", 
                                       round((max(peak_values())-mean(peak_values()))/sd_noise_estimate))),
           text.col = c("black", "green", "red", "blue"))
  })
  
  output$xcmsPlot <- renderPlot({
    chr <- Chromatogram(rtime = 1:length(peak_values()), intensity = peak_values())
    plot(findChromPeaks(chr, param = CentWaveParam()), main="")
    abline(v=xcms:::.getRtROI(intensity(chr), rtime(chr))[,3], col="red")
  })
}

# Run app ----
shinyApp(ui = ui, server = server)




