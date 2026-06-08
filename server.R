# server.R

server <- function(input, output, session) {
  addResourcePath("splicing_plots", file.path(getwd(), "R_FRASER"))
  
  # 1. Global Clear Trigger
  global_clear_trigger <- reactive({ input$clear_selections })
  
  # 2. Global IGV Ready Signal
  igv_ready_signal <- reactive({ input$igvReady })
  
  # 3. Launch Variants Module (FIXED: Added igv_ready argument)
  variants_Server("variants", 
                  clear_trigger = global_clear_trigger,
                  igv_ready = igv_ready_signal)
  
  # 4. Launch Patient Module
  patient_Server("patients", 
                 clear_trigger = global_clear_trigger,
                 igv_ready = igv_ready_signal)
}
