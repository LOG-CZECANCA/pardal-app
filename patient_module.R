# patient_module.R

patient_UI <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # Patient Lookup Search Box
      div(style = "margin-bottom: 20px;",
          h4("Patient Lookup", style = "margin-top: 0;"),
          selectizeInput(
            ns("selected_patient"),
            label = "Search / Select Patient ID(s):",
            choices = NULL, # Populated by the server for speed
            multiple = TRUE,
            options = list(placeholder = 'Type a patient ID...', maxItems = 5)
          )
      ),
      tags$hr(style = "border-color: #ddd; margin-top: 5px; margin-bottom: 15px;"),
      
      # Output View Buttons
      div(class = "output-radio-group",
          radioButtons(
            ns("patient_view_mode"),
            label = "Output",
            choices = c("Variants & Parallel Analysis" = "variants", 
                        "Allele Specific Expression" = "ase",
                        "Read Counts" = "counts",
                        "Expression Outliers (Outrider)" = "outrider",
                        "Splicing Outliers (Fraser)" = "fraser"), 
            selected = "variants"
          )
      ),
      
      tags$hr(),
      h4("Downloads"),
      
      # Transparent Download Buttons
      div(style = "display: flex; flex-direction: column; gap: 5px;",
          downloadButton(ns("dl_main"), "Variants (List)", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;"),
          downloadButton(ns("dl_parallel"), "Variants (Parallel)", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;"),
          downloadButton(ns("dl_ase"), "Allele Specific Expression", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;"),
          downloadButton(ns("dl_outrider"), "Outrider Outliers", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;"),
          downloadButton(ns("dl_fraser"), "Fraser Outliers", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;"),
          downloadButton(ns("dl_counts"), "Read Counts", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;"),
          tags$hr(style = "margin: 5px 0; border-color: #ddd;"),
          downloadButton(ns("dl_joined"), "Joined Variants Data", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;"),
          downloadButton(ns("dl_everything"), "Download ALL (ZIP)", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;")
      )
    ),
    
    mainPanel(
      width = 9,
      tags$head(tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min.js")),
      
      # --- A. VARIANTS VIEW ---
      conditionalPanel(
        condition = "input.patient_view_mode == 'variants'",
        ns = ns,
        do.call(tabsetPanel, c(id = ns("column_groups_patient"), make_group_tabs("patient", ns))),
        tags$hr(),
        
        div(style = "margin-bottom: 15px; background: #f9f9f9; padding: 15px; border-radius: 5px; border: 1px solid #ddd;",
            tags$strong("Active Filters:", style = "color: #555; display: block; margin-bottom: 10px;"),
            checkboxInput(ns("filter_gq"), "GQ=99 Filter", value = TRUE),
            checkboxInput(ns("filter_dp"), "DP>10 Filter", value = TRUE),
            checkboxInput(ns("filter_af"), "AF<0.01 Filter", value = TRUE),
            checkboxInput(ns("filter_clinsig"), "'Benign' Filter", value = TRUE),
            checkboxInput(ns("filter_feature"), "'Missing Feature' filter", value = TRUE)
        ),
        
        div(style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px;",
            h4("List of DNA/RNA Variants", style="margin: 0;"),
            actionButton(ns("btn_apply_global_filter"), "Apply Filters Globally", class="btn-sm", style="background-color: #FC65AB; color: white; border: none;")
        ),
        DTOutput(ns("patient_table_output")),
        tags$hr(),
        
        div(style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px;",
            h4("Parallel DNA/RNA Info", style="margin: 0;")
        ),
        DTOutput(ns("patient_parallel_table")),
        tags$hr(),
        
        div(
          style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px;",
          h4("IGV Viewer & Details (Variants)", style = "margin: 0;"),
          actionButton(inputId = ns("btn_screenshot_var"), label = "Save Image", icon = icon("camera"), class = "btn-primary btn-sm", onclick = "captureIGV('igv-capture-container-variants')")
        ),
        div(id = "igv-capture-container-variants", style = "background-color: white; padding: 10px; border-radius: 4px;", 
            div(id = "variant-details-container-var", class = "no-clear-click",
                wellPanel(style = "padding: 10px; background-color: #f9f9f9; border: 1px solid #ddd; margin-bottom: 10px;",
                          h5("Selected Event Information:", style = "margin-top:0px; font-weight:bold;"),
                          htmlOutput(ns("variant_details_box_var")))
            ),
            igvShinyOutput(ns("igvShiny_variants"), height = "600px")
        )
      ),
      
      # --- B. OTHER VIEWS ---
      conditionalPanel(
        condition = "input.patient_view_mode == 'ase'", ns = ns,
        div(style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
            h4("Allele Specific Expression (selected patient)", style="margin: 0;"),
            actionButton(ns("btn_apply_ase_filter"), "Apply Filters Globally", class="btn-sm", style="background-color: #FC65AB; color: white; border: none;")
        ),
        DTOutput(ns("patient_ase_table"))
      ),
      
      # --- COUNTS VIEW ---
      conditionalPanel(
        condition = "input.patient_view_mode == 'counts'", ns = ns,
        checkboxInput(ns("filter_low_counts"), "Filter out Read Counts <= 3", value = TRUE),
        div(style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
            h4("Gene Counts (selected patient)", style="margin: 0;"),
            actionButton(ns("btn_apply_counts_filter"), "Apply Filters Globally", class="btn-sm", style="background-color: #FC65AB; color: white; border: none;")
        ),
        DTOutput(ns("patient_counts_table")),
        
        tags$hr(style = "margin-top: 10px; margin-bottom: 10px;"),
        
        h4("GTEx Tissue Expression", style="margin-top: 10px;"), 
        div(class = "no-clear-click", style = "margin-bottom: 15px; border: 1px solid #eee; padding: 15px; border-radius: 4px; background-color: white; overflow-x: auto;", 
            plotly::plotlyOutput(ns("gtex_expression_plot"), height = "1000px")
        )
      ),

      conditionalPanel(
        condition = "input.patient_view_mode == 'outrider'", ns = ns,
        div(style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
            h4("Expression Outliers (Outrider) (selected patient)", style="margin: 0;"),
            actionButton(ns("btn_apply_outrider_filter"), "Apply Filters Globally", class="btn-sm", style="background-color: #FC65AB; color: white; border: none;")
        ),
        DTOutput(ns("patient_outrider_table"))
      ),
      
      conditionalPanel(
        condition = "input.patient_view_mode == 'fraser'", ns = ns,
        div(style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
            h4("Splicing Outliers (Fraser)", style="margin: 0;"),
            actionButton(ns("btn_apply_fraser_filter"), "Apply Filters Globally", class="btn-sm", style="background-color: #FC65AB; color: white; border: none;")
        ),
        div(class = "no-clear-click", verbatimTextOutput(ns("fraser_debug_info"))),
        DTOutput(ns("patient_fraser_table")),
        tags$hr(),
        div(style = "display: flex; align-items: center; background: #f5f5f5; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
            icon("magnifying-glass"), span("Zoom:", style="margin-left: 5px; font-weight: bold; margin-right: 15px;"),
            sliderInput(ns("plot_zoom"), label = NULL, min = 50, max = 250, value = 100, step = 10, width = "300px", post = "%")
        ),
        div(class = "no-clear-click", style = "overflow: auto; border: 1px solid #ddd; max-height: 800px; text-align: center; background: white; min-height: 100px;",
            uiOutput(ns("fraser_plot_output"))
        ),
        tags$hr(),
        div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px;",
            h4("IGV Viewer & Details (Splicing)", style = "margin: 0;"),
            actionButton(inputId = ns("btn_screenshot_fra"), label = "Save Image", icon = icon("camera"), class = "btn-primary btn-sm", onclick = "captureIGV('igv-capture-container-fraser')")
        ),
        div(id = "igv-capture-container-fraser", style = "background-color: white; padding: 10px; border-radius: 4px;", 
            div(id = "variant-details-container-fra", class = "no-clear-click",
                wellPanel(style = "padding: 10px; background-color: #f9f9f9; border: 1px solid #ddd; margin-bottom: 10px;",
                          h5("Selected Event Information:", style = "margin-top:0px; font-weight:bold;"),
                          htmlOutput(ns("variant_details_box_fra")))
            ),
            igvShinyOutput(ns("igvShiny_fraser"), height = "600px")
        )
      )
    )
  )
}

patient_Server <- function(id, clear_trigger = NULL, igv_ready = NULL) {
  moduleServer(id, function(input, output, session) {
    
    # Automatically populate the patient dropdown
    observe({
      if (exists("parallel_data") && !is.null(parallel_data)) {
        patients <- sort(unique(as.character(parallel_data$sample_id)))
        # Filter out empty or NA values just to be safe
        patients <- patients[!is.na(patients) & patients != ""]
        updateSelectizeInput(session, "selected_patient", choices = patients, server = TRUE)
      }
    })
    
    patient_table_proxy <- dataTableProxy('patient_table_output')
    patient_parallel_proxy <- dataTableProxy('patient_parallel_table')
    
    selected_patient_variant <- reactiveVal(NULL)
    selected_patient_parallel_variant <- reactiveVal(NULL)
    selected_patient_fraser_variant <- reactiveVal(NULL)
    
    selected_patient_ase_row <- reactiveVal(NULL)
    selected_patient_counts_row <- reactiveVal(NULL)
    selected_patient_outrider_row <- reactiveVal(NULL)
    selected_patient_outlier_bridge <- reactiveVal(NULL)
    
    selected_symbol_all <- reactiveVal(NULL)
    
    p_cf_main <- reactiveVal(NULL)
    p_cf_parallel <- reactiveVal(NULL)
    p_cf_ase <- reactiveVal(NULL)
    p_cf_counts <- reactiveVal(NULL)
    p_cf_outrider <- reactiveVal(NULL)
    p_cf_fraser <- reactiveVal(NULL)
    
    p_sync_main <- reactiveVal(NULL)
    p_sync_parallel <- reactiveVal(NULL)
    
    source("patient_reactives.R", local = TRUE)
    source("patient_observers.R", local = TRUE) 
    source("patient_outputs.R", local = TRUE) 
  })
}
