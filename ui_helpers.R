# ui_helpers.R

make_group_tabs <- function(id_prefix, ns = NULL) {
  if(is.null(ns)) ns <- function(x) x
  lapply(names(group_list), function(group_name) {
    group_id <- gsub(" ", "_", group_name)
    
    # We check which columns in this specific group are in our 'required' list
    default_selections <- intersect(group_list[[group_name]], required_columns)
    
    tabPanel(
      title = group_name,
      div(
        style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
        tags$h4(group_name),
        checkboxInput(inputId = ns(paste0("toggle_", id_prefix, "_", group_id)), label = NULL, value = FALSE)
      ),
      checkboxGroupInput(
        inputId = ns(paste0("group_", id_prefix, "_", group_id)),
        label = NULL,
        choices = group_list[[group_name]],
        selected = default_selections # Pre-selects the defaults
      )
    )
  })
}

app_css <- HTML("
/* Use box-shadow instead of background-color for DataTables v1.12+ */
table.dataTable tbody tr:not(.selected):hover,
table.dataTable tbody tr:not(.selected):hover td {
  box-shadow: inset 0 0 0 9999px #FD64AB !important;
  color: white !important;
}
table.dataTable tbody tr.selected,
table.dataTable tbody tr.selected td {
  box-shadow: inset 0 0 0 9999px #FD64AB !important;
  color: white !important;
}
table.dataTable.display tbody tr.selected:hover,
table.dataTable.display tbody tr.selected:hover > td,
table.dataTable.display tbody tr.selected:hover td {
  box-shadow: inset 0 0 0 9999px #FD64AB !important;
  color: white !important;
}

.output-radio-group .radio label {
  display: block;
  padding: 8px 12px;
  border-radius: 4px;
  background-color: transparent;
  color: black;
  border: 1px solid #ccc;
  cursor: pointer;
  transition: 0.2s;
  margin-bottom: 5px;
}

/* When selected: Pink Background, White Text */
.output-radio-group .radio label:has(input[type='radio']:checked) {
  background-color: #FC65AB;
  color: white;
  border-color: #FC65AB;
}

/* Hide the default circular radio dot */
.output-radio-group .radio input[type='radio'] {
  display: none;
}
")
# --- ROBUST CLICK LISTENER & IGV HANDLERS ---
app_js <- HTML("
    $(document).ready(function() {
      $(document).on('click', function(e) {
        
        if ($(e.target).closest('.no-clear-click').length > 0) {
            return;
        }

        if (
            $(e.target).closest('.dataTables_wrapper').length || 
            $(e.target).closest('input').length ||
            $(e.target).closest('select').length ||
            $(e.target).closest('.btn').length ||
            $(e.target).closest('.checkbox').length ||
            $(e.target).closest('.radio').length ||
            $(e.target).closest('.form-control').length ||
            $(e.target).closest('.selectize-control').length ||
            $(e.target).closest('.irs').length || 
            $(e.target).closest('.sidebar-panel').length || 
            $(e.target).closest('.nav-tabs').length ||
            $(e.target).closest('.igv-container').length ||
            $(e.target).closest('.igv-root-div').length ||
            $(e.target).closest('.html-widget').length ||
            $(e.target).closest('img').length 
        ) { 
            return; 
        }

        Shiny.setInputValue('clear_selections', Math.random());
      });
    });

    Shiny.addCustomMessageHandler('forceTableUpdate', function(tableId) {
      if (window[tableId + '_table']) {
        window[tableId + '_table'].ajax.reload(null, false);
      }
    });
    
    Shiny.addCustomMessageHandler('igv-nav', function(message) {
      var elementId = message.id; 
      var locus = message.locus;
      var container = document.getElementById(elementId);
      if (container && container.igvBrowser) {
        container.igvBrowser.search(locus);
      }
    });

    /* --- NEW: BULLETPROOF TRACK CLEARER --- */
    Shiny.addCustomMessageHandler('force-remove-igv-tracks', function(message) {
      var container = document.getElementById(message.id);
      if (container && container.igvBrowser) {
        message.trackNames.forEach(function(tname) {
          container.igvBrowser.removeTrackByName(tname);
        });
      }
    });
    
    window.captureIGV = function(containerId) {
      var targetId = containerId;
      if (!targetId) targetId = 'igv-capture-container-variants';
      
      var element = document.getElementById(targetId);
      
      if (!element) {
        console.error('Capture target not found: ' + targetId);
        return;
      }
      
      var nameDiv = element.querySelector('#igv-snapshot-filename');
      var fileName = 'IGV_Snapshot.png'; 
      
      if (nameDiv && nameDiv.innerText) {
        fileName = nameDiv.innerText.trim();
      }
      
      html2canvas(element, { 
        useCORS: true, 
        backgroundColor: '#ffffff'
      }).then(function(canvas) {
        var link = document.createElement('a');
        link.download = fileName;
        link.href = canvas.toDataURL();
        link.click();
      });
    };
")
