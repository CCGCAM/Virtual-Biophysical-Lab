
#rm(list= ls())

# 1. load the main libraries  -----


# 1. Install and load ToolsRTM from GitLab if missing -----------------------
if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  # Ensure remotes is available for GitLab installation
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_gitlab("caminoccg/toolsrtm", upgrade = "never")
}
library(ToolsRTM)
# Required packages

required_packages <- c("shiny",'expint','signal', 'textshaping',
                       "shinythemes", "ggplot2", "plotly", "reshape2", "tidyr","dplyr",'DT')

# 2.Check for missing packages and install them if necessary  -----

missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
# Load the libraries
lapply(required_packages, library, character.only = TRUE)

# 3.Define UI -----

ui <- navbarPage("Online reflectance simulator",theme = shinytheme("flatly"),

          tabPanel(title = "Virtual Biophysical Lab",
  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(width = 6,

                 fluidRow(h4("Leaf parameters"),

                          column(width = 4,

                                 sliderInput(inputId = "Cab",
                                             label = "Chlorophyll content",
                                             min = 0,
                                             max = 80,
                                             value = 40),
                                 sliderInput(inputId = "Cbrown",
                                             label = "Brown pigments",
                                             min = 0,
                                             max = 1,
                                             step = 0.05,
                                             value = 0.0),


                          ),

                          column(width = 4,
                                 sliderInput(inputId = "Car",
                                             label = "Carotenoid content",
                                             min = 0,
                                             max = 20,
                                             value = 8),
                                 sliderInput(inputId = "EWT",
                                             label = "Equivalent water thickness",
                                             min = 0.0001,
                                             max = 0.05,
                                             value = 0.009),


                          ),

                          column(width = 4,
                                 sliderInput(inputId = "Anth",
                                             label = "Anthocyanin Content",
                                             min = 0,
                                             max = 7,
                                             step = 0.2,
                                             value = 2),
                                 sliderInput(inputId = "CBC",
                                             label = "Carbon-based constituent",
                                             min =  0.0001,
                                             max = 0.03,
                                             value = 0.005),

                          ),
                          column(width = 4,
                                 sliderInput(inputId = "N",
                                             label = "Structure parameter N",
                                             min = 1,
                                             max = 3,
                                             step = 0.2,
                                             value = 1.5),


                          ),


                          column(width = 4,

                                 sliderInput(inputId = "Prot",
                                             label = "Leaf protein content",
                                             min = 0.0001,
                                             max = 0.03,
                                             value = 0.0045)
                          )),

                 fluidRow(h4("Canopy parameters"),

                          column(width = 4,

                                 sliderInput(inputId = "LAI",
                                             label = "Leaf Area Index",
                                             min = 0.001,
                                             max = 10,
                                             step = 0.1,
                                             value = 4)),

                          column(width = 4,

                                 sliderInput(inputId = "hspot",
                                             label = "Hotspot parameter",
                                             min = 0,
                                             max = 1,
                                             value = 0.01)),

                          column(width = 4,

                                 selectInput(inputId = "TypeLIDF",
                                             label = "Type LIDF",
                                             choices = list(Planophile = "plano",
                                                            Erectophile = "erecto",
                                                            Plagiophile = "plagio",
                                                            Extremophile = "extremo",
                                                            Uniform = 'uniform',
                                                            Spherical = "sph"),
                                             selected = "sph"))
                 ),

                 fluidRow(h4("Soil parameters"),

                          column(width = 4,

                                 sliderInput(inputId = "psoil",
                                             label = "Soil brightness",
                                             min = 0,
                                             max = 1,
                                             value = 0.5)),
                          # Checkbox input for including soil reflectance
                          column(width = 4,
                                 checkboxInput(inputId = "include_soil_reflectance",
                                               label = "Include soil reflectance",
                                               value = F)  # Default to TRUE or FALSE as per your preference
                          )
                 ),

                 fluidRow(h4("Angle parameters"),

                          column(width = 4,

                                 sliderInput(inputId = "tts",
                                             label = "Solar zenith angle",
                                             min = 0,
                                             max = 90,
                                             value = 30)),
                          column(width = 4,

                                 sliderInput(inputId = "tto",
                                             label = "Observer zenith angle",
                                             min = 0,
                                             max = 90,
                                             value = 10)),

                          column(width = 4,

                                 sliderInput(inputId = "psi",
                                             label = "Relative azimuth angle",
                                             min = 0,
                                             max = 180,
                                             value = 0))
                          ),

                 # Add options to visualize reflectance and transmittance at leaf level
                 fluidRow(h4("Leaf Optical Properties"),

                          column(width = 4,
                                 checkboxInput(inputId = "show_reflectance",
                                               label = "Show Reflectance at Leaf Level",
                                               value = FALSE)),


                          column(width = 4,
                                 checkboxInput(inputId = "show_transmittance",
                                               label = "Show Transmittance at Leaf Level",
                                               value = FALSE)),
                          column(width = 4,
                                 checkboxInput(inputId = "show_diff",
                                               label = "Show spectral differences",
                                               value = FALSE)),

                 ),     fluidRow(
                   column(width = 6,  # Adjust the width as needed
                          div(style = "overflow-y: auto; height: 350px;",
                              downloadButton("downloadData", "Download PROSAIL Data"))
                   ),
                   #column(width = 6,  # Adjust the width as needed
                          #div(style = "overflow-y: auto; height: 350px;",
                           #  checkboxInput(inputId = "accum_spectra",
                            #                label = "Accumulative spectra for downloading",
                            #                value = FALSE))
                   #)


                 )
    ),

    # Main panel for displaying outputs ----
    mainPanel(width = 6,


              # Create a tabset panel
              tabsetPanel(id = "tabs",
                tabPanel("Simulations",

                         h4("Simulated Spectrum"),
                         div(style = "overflow-y: auto; height: 400px;",  # Set a fixed height and enable vertical scrolling
                         #    plotOutput(outputId = "simulations") ),
                         # UI ----
                         plotlyOutput("simulations", height = "420px") ),
                         div(style = "overflow-y: auto; height: 400px;",  # Combined plot (reflectance and transmittance)
                             plotOutput(outputId = "refl_trans_plot")   )
                ),
                tabPanel("Spectral Differences",
                         div(style = "overflow-y: auto; height: 800px;",  # Set a fixed height for the diff plot
                             plotOutput(outputId = "diff_plot",height = "600px"))
                ),
                tabPanel("Absorption coefficients",
                         div(style = "overflow-y: auto; height: 350px;",  # Set a fixed height for the diff plot
                             plotOutput(outputId = "abs_plot1",height = "350px")),
                         div(style = "overflow-y: auto; height: 350px;",  # Set a fixed height for the diff plot
                             plotOutput(outputId = "abs_plot2",height = "350px"))
                )

              )

    )
  )
)
)
# 4. Define server logic -----

server <- function(input, output,session) {


  # Define the reactive expression for the leaf parameters ---------------------------------------------------
  parameters_model <- reactive({
    params <- c(
      N = as.numeric(input$N),
      Cab = as.numeric(input$Cab),Car = as.numeric(input$Car),Anth = as.numeric(input$Anth),
      Cbrown = as.numeric(input$Cbrown),EWT = as.numeric(input$EWT),
      LMA = 0,                # Assuming LMA is constant, no need for conversion
      alpha = 40,             # Same with alpha
      Prot = as.numeric(input$Prot),CBC = as.numeric(input$CBC),
      LAI = as.numeric(input$LAI),
      TypeLidf = 1,           # Assuming 1 for considering LIDFb=0. and LIDFa in degrees.
      hspot = as.numeric(input$hspot),
      tts = as.numeric(input$tts),tto = as.numeric(input$tto),psi = as.numeric(input$psi),
      psoil = as.numeric(input$psoil))

    return(params)
  })


  # Observe any slider input change and switch to 'Simulations' tab
  observeEvent({
    # List all input values that should trigger the tab switch
    list(input$N, input$Cab, input$Car, input$Cbrown, input$Prot, input$CBC, input$EWT,
         input$LAI, input$hspot,
         input$TypeLIDF,
         input$tts, input$tto, input$psi, input$psoil)
  }, {
    updateTabsetPanel(session, "tabs", selected = "Simulations")
  })


  # Reactive expression for the LUT data ---------------------------------------------------
  lut_data.sim <- reactive({
    # Get the selected leaf and canopy parameters from the reactive `params` expression
    lut.to_sim <- data.frame(t(as.data.frame(parameters_model())))
    row.names(lut.to_sim) <- NULL


    # Call the Simulations
    data <- ToolsRTM::dataSpec_PDB
    Rsoil.dry  <- data[,11]  # rsoil1 = dry soil
    Rsoil.wet <- data[,12]   # rsoil2 = wet soil
    print(lut.to_sim$psoil)
    psoil <- lut.to_sim[1,'psoil']

    # soil factor (psoil=0: wet soil / psoil=1: dry soil)
    rsoil <- psoil * Rsoil.dry + (1 - psoil) * Rsoil.wet

    return(list(lut.to_sim, rsoil))
  })


  # Observe the checkbox input to switch to the "Spectral Differences" tab
  observeEvent(input$show_diff, {
    if (input$show_diff) {
      updateTabsetPanel(session, "tabs", selected = "Spectral Differences")
    }
  })

  # Observe the checkbox input to switch to the "Spectral Differences" tab
  observeEvent(input$show_reflectance, {
    if (input$show_reflectance) {
      updateTabsetPanel(session, "tabs", selected = "Simulations")
    }
  })

  # Observe the checkbox input to switch to the "Spectral Differences" tab
  observeEvent(input$show_transmittance, {
    if (input$show_transmittance) {
      updateTabsetPanel(session, "tabs", selected = "Simulations")
    }
  })

  # Reactive expression for reflectance data ---------------------------------------------------
  reflectance_data <- reactive({


    df.LUT <- lut_data.sim()[[1]]
    row.names(df.LUT) <- NULL
    print(df.LUT)

    rsoil <- lut_data.sim()[[2]]

    # Update LIDF parameters based on input
    df.LUT$LIDFa <- switch(input$TypeLIDF,
                           plano = 1,
                           erecto = -1,
                           plagio = 0,
                           extremo = 0,
                           uniform = 0,
                           sph = -0.35)

    df.LUT$LIDFb <- switch(input$TypeLIDF,
                           plano = 0,
                           erecto = 0,
                           plagio = -1,
                           extremo = 1,
                           uniform = 0,
                           sph = -0.15)


    print(df.LUT)


    # Compute PROSPECT-PRO simulation as leaf model ---------------------------------------------------
    sim_leaf_values <- ToolsRTM::prospect_PRO(N=df.LUT$N,Cab=df.LUT$Cab, Car=df.LUT$Car, Anth=df.LUT$Anth,Cbrown=df.LUT$Cbrown,
                                            EWT=df.LUT$EWT,LMA=df.LUT$LMA,alpha=40,
                                            Prot=df.LUT$Prot,CBC=df.LUT$CBC)
    wavelength <- sim_leaf_values[[1]]
    rho	 <- 	sim_leaf_values[[2]]
    tau	 <- 	sim_leaf_values[[3]]

    # Compute foursail simulation adding PROSAIL as leaf model  ---------------------------------------------------
    reflectance_values <- ToolsRTM::foursail(inputLUT=df.LUT[1,],rsoil=rsoil, LeafModel = 'PROSPECT-PRO')
    rdot<-reflectance_values[[1]]
    rsot<-reflectance_values[[2]]
    showNotification("Forward simulation done successfully.", type = "message")


    # Compute the Bidirectional Reflectance Factor (BRDF) ---------------------------------------------------
    reflectance_values<- ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=df.LUT[1,'tts'],data.light=ToolsRTM::dataSpec_PDB)
    #showNotification("Forward simulation with BRDF model successfully applied.", type = "message")

    reflectance_df <- data.frame(wavelength = wavelength, reflectance = reflectance_values)


    # Plot reflectance/ transmittance at leaf level  ---------------------------------------------------
    output$refl_trans_plot <- renderPlot({
      # Check if either reflectance or transmittance is selected
      if (input$show_reflectance || input$show_transmittance) {

        # Initialize an empty dataframe for combined plotting
        df_combined <- data.frame()

        # Add reflectance data if selected
        if (input$show_reflectance) {
          df_combined <- rbind(
            df_combined,
            data.frame(wavelength = wavelength, value = rho, variable = "Reflectance")
          )
        }

        # Add transmittance data if selected
        if (input$show_transmittance) {
          df_combined <- rbind(
            df_combined,
            data.frame(wavelength = wavelength, value = tau, variable = "Transmittance")
          )
        }
        # Show a notification when the simulation data is processed
        showNotification("Leaf-level simulation successfully completed.", type = "message")


        # Plot both reflectance and transmittance in a single plot
        ggplot(df_combined, aes(x = wavelength, y = value, color = variable)) +
          geom_line(linewidth = 1) +
          facet_grid(variable ~ ., scales = "free_y") +  # Create facets for reflectance and transmittance
          scale_color_manual(values = c("Reflectance" = "forestgreen", "Transmittance" = "navyblue")) +
          labs(x = "Wavelength (nm)", y = NULL) +  # No common y-axis label
          theme_bw() +
          theme(
            legend.position = 'none',
            strip.text = element_text(face = "bold", size = 14),  # Bold facet labels
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(face = "bold", size = 14),
            axis.text.y = element_text(hjust = 0.5, size = 12, face = "bold"),
            axis.text.x = element_text(hjust = 0.5, size = 12, face = "bold"),
            panel.grid.major = element_blank(),  # Optional: Remove grid lines
            panel.grid.minor = element_blank())
      }
    })


    # Plot spectral differences between canopy and leaf  ---------------------------------------------------

    output$diff_plot <- renderPlot({
      if (input$show_diff) {  # Add checkbox input for showing differences

        # Create the dataframe for differences
        df.diff<- data.frame(
          wavelength = wavelength,
          refl.canopy = reflectance_values,
          refl.leaf = rho,
          diff = abs(reflectance_values - rho),
          diff_sqrt = (reflectance_values - rho)^2,
          percent_diff = (reflectance_values - rho) / rho * 100)

        # Reshape the data for plotting with facets
        df.melted <- reshape2::melt(df.diff, id.vars = "wavelength",
                                    measure.vars = c("diff", 'diff_sqrt',"percent_diff"),
                                    variable.name = "variable",
                                    value.name = "value")


        # Plot the difference
        plot.dif <-ggplot(df.diff, aes(x = wavelength, y = diff)) +
          geom_line(color = "purple", linewidth = 1) +
          labs(x = "Wavelength (nm)", y = "Difference (Canopy - Leaf Reflectance)") +
          theme_bw() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(face = "bold", size = 14),
            axis.text.y = element_text(hjust = 0.5, size = 12, face = "bold"),
            axis.text.x = element_text(hjust = 0.5, size = 12, face = "bold"))



        # Define custom labels for facets
        custom_labels <- c(
          "diff" = "Absolute Difference",
          'diff_sqrt' = 'Squared Difference',
          "percent_diff" = "Difference in %"
        )

        # Plot using facets
        plot.dif <- ggplot(df.melted, aes(x = wavelength, y = value, color = variable)) +
          geom_line(linewidth = 1) +
          facet_grid(variable ~ ., scales = "free_y", labeller = as_labeller(custom_labels)) +  # Create facets with custom labels
          scale_color_manual(values = c("diff" = "gray6", "percent_diff" = "dodgerblue4", 'diff_sqrt' = 'darkolivegreen4')) +
          labs(x = "Wavelength (nm)",
               y = NULL,
               title = "Difference Between Canopy and Leaf Reflectance") +
          theme_bw() +
          theme(
            legend.position = 'none',
            strip.text = element_text(face = "bold", size = 14),  # Bold facet labels
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(face = "bold", size = 14),
            axis.text.y = element_text(hjust = 0.5, size = 12, face = "bold"),
            axis.text.x = element_text(hjust = 0.5, size = 12, face = "bold"),
            panel.grid.major = element_blank(),  # Optional: Remove grid lines
            panel.grid.minor = element_blank()
          )


        # Show notification when the plot is generated
        showNotification("Difference between canopy and leaf reflectance successfully generated.", type = "message")
        #print(plot.dif)
        print(plot.dif)
      }
    })


    return(reflectance_df)



  })




  # Plot PROSAIL at canopy level  ---------------------------------------------------
  output$simulations <- plotly::renderPlotly({

    df.plot <- reflectance_data()

    # Base ggplot
    plot_spectra <- ggplot(df.plot, aes(x = wavelength, y = reflectance)) +
      # Highlight regions with semi-transparent rectangles
      annotate("rect", xmin = 400, xmax = 700, ymin = -0.05, ymax = 0.7,
               fill = "skyblue", alpha = 0.1) +   # VIS
      annotate("rect", xmin = 700, xmax = 1300, ymin = -0.05, ymax = 0.7,
               fill = "lightgreen", alpha = 0.1) + # NIR
      annotate("rect", xmin = 1300, xmax = 2500, ymin = -0.05, ymax = 0.7,
               fill = "orange", alpha = 0.1) +   # SWIR
      geom_line(aes(color = "RTM"), linewidth = 1) +
      labs(x = "Wavelength (nm)", y = "Reflectance", color = NULL) +
      # Add region labels
      annotate("text", x = 550,  y = 0.7, label = "VIS",  size = 5, fontface = "bold", color = "blue") +
      annotate("text", x = 1000, y = 0.7, label = "NIR",  size = 5, fontface = "bold", color = "darkgreen") +
      annotate("text", x = 1900, y = 0.7, label = "SWIR", size = 5, fontface = "bold", color = "darkorange") +

      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title  = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold")
      ) +
      scale_color_manual(values = c("RTM" = "black"))

    # Soil (opcional)
    if (isTRUE(input$include_soil_reflectance)) {
      rsoil <- lut_data.sim()[[2]]
      df.soil <- data.frame(wavelength = df.plot$wavelength, reflectance = rsoil)

      plot_spectra <- plot_spectra +
        geom_line(data = df.soil, aes(x = wavelength, y = reflectance, color = "Soil"), linewidth = 1) +
        scale_color_manual(values = c("RTM" = "black", "Soil" = "goldenrod4"))

      showNotification("Soil reflectance spectrum is added successfully.", type = "message")
    }


    # Convert to Plotly and remove legend
    p <- plotly::ggplotly(plot_spectra, tooltip = c("x", "y")) |>
      plotly::layout(showlegend = FALSE) |>
      plotly::config(displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d"))

    p
  })





  output$abs_plot1 <- renderPlot({

    df.plot <-ToolsRTM::dataSpec_PRO
    head(df.plot)
    # Subset the data for the wavelength range 400-900 nm
    df.plot <- df.plot[df.plot$wave >= 400 & df.plot$wave <= 900, ]

    # Reshape df.plot to long format for easier plotting
    df_long <- melt(df.plot, id.vars = "wave",
                    measure.vars = c("spA_Cab", "spA_Car", "spA_Atn", "spA_Br"),
                    variable.name = "Trait",
                    value.name = "SpecificLeafArea")
    df_long$Trait <- factor(df_long$Trait,
                            levels = c("spA_Cab", "spA_Car", "spA_Atn", "spA_Br"),
                            labels = c("Chlorophyll content", "Carotenoids content", "Anthocyanin content", "Brown pigments"))



    plot_spectra <- ggplot(df_long, aes(x = wave, y = SpecificLeafArea, color = Trait)) +
      geom_line(linewidth=1) +
      labs(title = "",
           x = "Wavelength (nm)",
           y = "Specific Absorption") +
      theme_bw() +
      theme(legend.position = 'right',
            legend.box.background = element_rect(color = "black",linewidth=1),
            plot.title = element_text(hjust = 0.5, size=14,face="bold"),
            axis.title = element_text(face="bold", size=14),
            legend.text = element_text(face="bold", size=12),
            axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
            axis.text.x=element_text(hjust = 0.5, size=12,face="bold"),
            legend.title=element_blank()) +
      scale_color_manual(values = c("Chlorophyll content" = "darkolivegreen",
                                    "Carotenoids content" = "darkorange3",
                                    "Anthocyanin content" = "firebrick3",
                                    "Brown pigments" = "gold4"))
    #  scale_color_brewer(palette = "Set1") # Optional: Choose a color palette

    print(plot_spectra)



    showNotification("Absorption spectrum done successfully.", type = "message")

  })




  output$abs_plot2 <- renderPlot({

    df.plot <- ToolsRTM::dataSpec_PRO
    # Reshape df.plot to long format for easier plotting
    df_long <- melt(df.plot, id.vars = "wave",
                     measure.vars = c( "spA_Cw", "spA_Protein", "spA_nonPro"), # "spA_Cm"
                     variable.name = "Trait",
                     value.name = "SpecificLeafArea")

    # Change legend labels
    df_long$Trait <- factor(df_long$Trait,
                            levels = c("spA_Cw", "spA_Protein", "spA_nonPro"), # "spA_Cm"
                            labels = c("Water content",  "Protein content", "Carbon-based compounds")) #"Dry matter content",


    plot_spectra <- ggplot(df_long, aes(x = wave, y = SpecificLeafArea, color = Trait)) +
      geom_line(linewidth=1) +
      labs(title = "",
           x = "Wavelength (nm)",
           y = "Specific Absorption") +
      theme_bw() +
      theme(legend.position = 'right',
            legend.box.background = element_rect(color = "black",linewidth=1),
            plot.title = element_text(hjust = 0.5, size=14,face="bold"),
            axis.title = element_text(face="bold", size=14),
            legend.text = element_text(face="bold", size=12),
            axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
            axis.text.x=element_text(hjust = 0.5, size=12,face="bold"),
            legend.title=element_blank()) +
      scale_color_manual(values = c("Water content" = "dodgerblue4",
                                  #  "Dry matter content" = "goldenrod4",
                                    "Protein content" = "firebrick3",
                                    "Carbon-based compounds" = "khaki3"))

    print(plot_spectra)


  })




  # Reactive value to track download success
  download_success <- reactiveVal(FALSE)
  # Create a download link
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("prosail_spectrum_", Sys.Date(), ".csv", sep = "")  # Added date to filename for uniqueness
    },
    content = function(file) {
      # Get the reflectance data
      df <- reflectance_data()  # Assuming this function retrieves the current reflectance data
      # Transpose the data
      df_transposed <- df %>%
        pivot_wider(names_from = wavelength, values_from = reflectance)
      # Write the accumulated data to a CSV file
      write.csv(reflectance_data(), file, row.names = FALSE)

    })
  # Set download success to TRUE
  download_success(TRUE)


  # Observe the download success reactive value to show notification
  observeEvent(download_success(), {
    if (download_success()) {
      showNotification("Downloaded PROSAIL data", type = "message")
      # Reset the reactive value to FALSE after showing the notification
      download_success(FALSE)
    }
  })

}

# 5. Create the Shiny app object -----

shinyApp(ui = ui, server = server)

