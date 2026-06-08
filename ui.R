# ui.R
ui <- fluidPage(
  tags$head(
    # 1. CSS for the Full-Screen Loading Overlay
    tags$style(HTML("
      #loading-screen {
        position: fixed;
        top: 0; left: 0; width: 100%; height: 100%;
        background-color: #ffffff;
        z-index: 9999;
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
      }
      .loader-spinner {
        border: 8px solid #f3f3f3;
        border-top: 8px solid #FC65AB; /* Matches your pink theme */
        border-radius: 50%;
        width: 70px;
        height: 70px;
        animation: spin 1s linear infinite;
      }
      @keyframes spin {
        0% { transform: rotate(0deg); }
        100% { transform: rotate(360deg); }
      }
    ")),
    
    # 2. JavaScript to hide the loading screen when the app is completely ready
    tags$script(HTML("
      /* 'shiny:idle' fires exactly when Shiny finishes all initial calculations and rendering */
      $(document).on('shiny:idle', function() {
        $('#loading-screen').fadeOut(500); /* Fades out over half a second */
      });
    "))
  ),
  
  # 3. The Loading Screen HTML
  div(id = "loading-screen",
      div(class = "loader-spinner"),
      h3("Loading Genomic Data...", style = "color: #555; margin-top: 20px; font-family: sans-serif;"),
      h5("Please wait while tables are assembled", style = "color: #888;")
  ),
  
  # --- Your existing UI elements below ---
  tags$style(app_css),
  tags$script(app_js),
  
  div(
    style = "text-align: center; padding: 0px;",
    img(src = "pardal.jpg", height = "150px", alt = "Header Image")
  ),
  
  tabsetPanel(
    id = "view_mode",
    
    tabPanel(
      "All Variants",
      variants_UI("variants")
    ),
    
    tabPanel(
      "Patient View",
      patient_UI("patients")
    )
  )
)
