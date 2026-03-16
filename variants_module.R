# variants_module.R

variants_UI <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # 🌟 Wrapped in our new CSS class to make them look like buttons
      div(class = "output-radio-group",
          radioButtons(
            ns("all_variants_view"),
            label = "Output",
            choices = c("Variants & Parallel Analysis" = "variants", 
                        "Allele Specific Expression" = "ase",
                        "Read Counts" = "counts",
                        "Expression Outliers (Outrider)" = "outrider",
                        "Splicing Outliers (Fraser)" = "fraser",
                        "Analyzed Sample Pairs" = "samples_list",       
                        "Sample Identity" = "sample_identity"), 
            selected = "variants"
          )
      ),
      tags$hr(),
      h4("Downloads"),
      
      # 🌟 Download Buttons changed to transparent background & black text
      div(style = "display: flex; flex-direction: column; gap: 5px;",
          downloadButton(ns("dl_all_main"), "Variants (Main)", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;"),
          downloadButton(ns("dl_all_parallel"), "Variants (Parallel)", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;"),
          downloadButton(ns("dl_ase"), "Allele Specific Expression", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;"),
          downloadButton(ns("dl_outrider"), "Outrider Outliers", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;"),
          downloadButton(ns("dl_fraser"), "Fraser Outliers", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;"),
          downloadButton(ns("dl_counts"), "Read Counts", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;"),
          tags$hr(style="margin: 5px 0; border-color: #ddd;"),
          downloadButton(ns("dl_all_joined"), "Joined Variants Data", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;"),
          downloadButton(ns("dl_everything"), "Download ALL (ZIP)", class = "btn-sm", style = "background-color: transparent; color: black; border: 1px solid #ccc;")
      )
    ),
    
    mainPanel(
      width = 9,
      tags$head(tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min.js")),
      
      conditionalPanel(
        condition = "input.all_variants_view == 'variants'",
        ns = ns,
        do.call(tabsetPanel, c(id = ns("column_groups"), make_group_tabs("all", ns))),
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
        DTOutput(ns("table_output_highlighted")),
        tags$hr(),
        div(style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px;",
            h4("Parallel DNA/RNA Info", style="margin: 0;")
        ),
        DTOutput(ns("parallel_table")),
        tags$hr(),
        div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px;",
            h4("IGV Viewer & Details", style = "margin: 0;"),
            actionButton(ns("btn_screenshot_var_global"), "Save Image", icon = icon("camera"), class = "btn-primary btn-sm", onclick = "captureIGV('igv-capture-container-variants-global')")
        ),
        div(id = "igv-capture-container-variants-global", style = "background-color: white; padding: 10px; border-radius: 4px;", 
            div(id = "variant-details-container-var-global", class = "no-clear-click",
                wellPanel(style = "padding: 10px; background-color: #f9f9f9; border: 1px solid #ddd; margin-bottom: 10px;",
                          h5("Selected Event Information:", style = "margin-top:0px; font-weight:bold;"),
                          htmlOutput(ns("variant_details_box_var")))
            ),
            igvShinyOutput(ns("igvShiny_variants_global"), height = "600px")
        )
      ),
      
      conditionalPanel(condition = "input.all_variants_view == 'ase'", ns = ns, div(style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;", h4("Allele Specific Expression (all patients)"), actionButton(ns("btn_apply_ase_filter"), "Apply Filters Globally", class="btn-sm", style="background-color: #FC65AB; color: white; border: none;")), DTOutput(ns("ase_table"))),
      
      # 🌟 CHECKBOX MOVED HERE TO BE VISIBLE IN COUNTS VIEW
      conditionalPanel(condition = "input.all_variants_view == 'counts'", ns = ns, 
                       checkboxInput(ns("filter_low_counts"), "Filter out Read Counts <= 3", value = TRUE),
                       div(style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;", h4("Read Counts"), actionButton(ns("btn_apply_counts_filter"), "Apply Filters Globally", class="btn-sm", style="background-color: #FC65AB; color: white; border: none;")), DTOutput(ns("counts_table")), tags$hr(), h4("GTEx Tissue Expression"), div(class = "no-clear-click", style = "margin-bottom: 15px; border: 1px solid #eee; padding: 15px; border-radius: 4px; background-color: white; overflow-x: auto;", plotly::plotlyOutput(ns("gtex_expression_plot"), height = "1000px"))),
      
      conditionalPanel(condition = "input.all_variants_view == 'outrider'", ns = ns, div(style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;", h4("Expression Outliers (Outrider)"), actionButton(ns("btn_apply_outrider_filter"), "Apply Filters Globally", class="btn-sm", style="background-color: #FC65AB; color: white; border: none;")), DTOutput(ns("outrider_table"))),
      conditionalPanel(condition = "input.all_variants_view == 'fraser'", ns = ns, div(style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;", h4("Splicing Outliers (Fraser)"), actionButton(ns("btn_apply_fraser_filter"), "Apply Filters Globally", class="btn-sm", style="background-color: #FC65AB; color: white; border: none;")), verbatimTextOutput(ns("fraser_debug_info")), DTOutput(ns("fraser_table"))),
      conditionalPanel(condition = "input.all_variants_view == 'samples_list'", ns = ns, h4("Analyzed Sample Pairs"), DTOutput(ns("samples_pairs_table"))),
      conditionalPanel(condition = "input.all_variants_view == 'sample_identity'", ns = ns, h4("Sample Identity"), DTOutput(ns("sample_identity_table")))
    )
  )
}

variants_Server <- function(id, clear_trigger = NULL, igv_ready = NULL) {
  moduleServer(id, function(input, output, session) {
    selected_variant_main <- reactiveVal(NULL)
    selected_variant_parallel <- reactiveVal(NULL)
    selected_fraser_variant_all <- reactiveVal(NULL)
    
    selected_ase_row <- reactiveVal(NULL)
    selected_counts_row <- reactiveVal(NULL)
    selected_outrider_row <- reactiveVal(NULL)
    selected_outlier_bridge <- reactiveVal(NULL) 
    
    selected_gene_all <- reactiveVal(NULL) 
    selected_sample_all <- reactiveVal(NULL)
    selected_symbol_all <- reactiveVal(NULL)
    
    v_cf_main <- reactiveVal(NULL)
    v_cf_parallel <- reactiveVal(NULL)
    v_cf_ase <- reactiveVal(NULL)
    v_cf_counts <- reactiveVal(NULL)
    v_cf_outrider <- reactiveVal(NULL)
    v_cf_fraser <- reactiveVal(NULL)
    
    v_sync_main <- reactiveVal(NULL)
    v_sync_parallel <- reactiveVal(NULL)
    
    source("variants_reactives.R", local = TRUE)
    source("variants_observers.R", local = TRUE)
    source("variants_outputs.R", local = TRUE)
  })
}