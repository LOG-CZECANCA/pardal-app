# patient_observers.R

if (!exists("loaded_igv_tracks_var", inherits = FALSE)) {
  loaded_igv_tracks_var <- reactiveVal(character(0))
  loaded_igv_tracks_fra <- reactiveVal(character(0))
  last_selected_var <- reactiveVal(NULL)
  last_selected_fra <- reactiveVal(NULL)
  current_igv_state_var <- reactiveVal(NULL)
  current_igv_state_fra <- reactiveVal(NULL)
}

is_actively_searched <- function(g, c) {
  has_g <- !is.null(g) && length(g) > 0 && g != ""
  has_c <- !is.null(c) && length(c) > 0 && any(c != "")
  return(has_g || has_c)
}

observe({
  current_selections <- as.character(unlist(lapply(names(group_list), function(g) { 
    input[[paste0("group_patient_", gsub(" ", "_", g))]] 
  })))
  
  old_order <- isolate(user_col_order_patient())
  kept_cols <- old_order[old_order %in% current_selections]
  new_cols <- setdiff(current_selections, old_order)
  
  user_col_order_patient(c(kept_cols, new_cols))
})

lapply(names(group_list), function(group_name) {
  group_id <- gsub(" ", "_", group_name)
  observeEvent(input[[paste0("toggle_patient_", group_id)]], {
    new_selections <- if (isTRUE(input[[paste0("toggle_patient_", group_id)]])) {
      unname(as.character(group_list[[group_name]]))
    } else { character(0) }
    updateCheckboxGroupInput(session, inputId = paste0("group_patient_", group_id), selected = new_selections)
  })
})

# --- 1. LOCAL AUTO-SYNC (Protected with identical() to prevent Defiltering Glitch) ---
observeEvent(input$patient_table_output_rows_all, {
  req(input$patient_table_output_rows_all)
  if (is_actively_searched(input$patient_table_output_search, input$patient_table_output_search_columns)) {
    df <- isolate(selected_data_patient())
    idx <- input$patient_table_output_rows_all
    idx <- idx[!is.na(idx) & idx > 0 & idx <= nrow(df)] 
    if (length(idx) > 0) {
      new_sync <- df[idx, , drop=FALSE]
      if (!identical(new_sync, p_sync_main())) p_sync_main(new_sync)
    } else { if (!is.null(p_sync_main())) p_sync_main(NULL) }
  } else { if (!is.null(p_sync_main())) p_sync_main(NULL) }
}, ignoreNULL = FALSE, ignoreInit = TRUE)

observeEvent(input$patient_parallel_table_rows_all, {
  req(input$patient_parallel_table_rows_all)
  if (is_actively_searched(input$patient_parallel_table_search, input$patient_parallel_table_search_columns)) {
    df <- isolate(patient_filtered_parallel_data())
    idx <- input$patient_parallel_table_rows_all
    idx <- idx[!is.na(idx) & idx > 0 & idx <= nrow(df)] 
    if (length(idx) > 0) {
      new_sync <- df[idx, , drop=FALSE]
      if (!identical(new_sync, p_sync_parallel())) p_sync_parallel(new_sync)
    } else { if (!is.null(p_sync_parallel())) p_sync_parallel(NULL) }
  } else { if (!is.null(p_sync_parallel())) p_sync_parallel(NULL) }
}, ignoreNULL = FALSE, ignoreInit = TRUE)


# --- 2. GLOBAL SYNC BUTTONS ---
broadcast_global_filter <- function(src, data) {
  payload <- list(src = src, data = data)
  
  if (src == "variants") {
    p_cf_main(NULL)
    p_cf_parallel(NULL)
  } else {
    if(src != "main") p_cf_main(payload)
    if(src != "parallel") p_cf_parallel(payload)
  }
  
  p_cf_ase(if (src == "ase") NULL else payload)
  p_cf_counts(if (src == "counts") NULL else payload)
  p_cf_outrider(if (src == "outrider") NULL else payload)
  p_cf_fraser(if (src == "fraser") NULL else payload)
}

observeEvent(input$btn_apply_global_filter, {
  df_main <- isolate(selected_data_patient())
  idx_main <- input$patient_table_output_rows_all
  if (is.null(idx_main)) idx_main <- seq_len(nrow(df_main))
  idx_main <- idx_main[!is.na(idx_main) & idx_main > 0 & idx_main <= nrow(df_main)]
  visible_main <- df_main[idx_main, , drop=FALSE]
  
  df_par <- isolate(patient_filtered_parallel_data())
  idx_par <- input$patient_parallel_table_rows_all
  if (is.null(idx_par)) idx_par <- seq_len(nrow(df_par))
  idx_par <- idx_par[!is.na(idx_par) & idx_par > 0 & idx_par <= nrow(df_par)]
  visible_par <- df_par[idx_par, , drop=FALSE]
  
  if(nrow(visible_main) > 0 && nrow(visible_par) > 0) {
    matching_par <- dplyr::semi_join(visible_par, visible_main, by = intersect(colnames(visible_par), c("chr", "pos", "ref", "alt", "strand")))
    if(nrow(matching_par) > 0) {
      broadcast_global_filter("variants", matching_par)
      showNotification("Global filters safely applied. Variants tables kept intact.", type="message")
    } else { showNotification("No overlapping variants found between tables.", type="warning") }
  } else { showNotification("No variants visible to apply.", type="warning") }
})

observeEvent(input$btn_apply_ase_filter, { df <- patient_filtered_ase(); visible <- df[input$patient_ase_table_rows_all, , drop=FALSE]; if(nrow(visible) > 0) { broadcast_global_filter("ase", visible); showNotification("Global filters applied from ASE.", type="message") } })
observeEvent(input$btn_apply_counts_filter, { df <- patient_filtered_counts(); visible <- df[input$patient_counts_table_rows_all, , drop=FALSE]; if(nrow(visible) > 0) { broadcast_global_filter("counts", visible); showNotification("Global filters applied from Gene Counts.", type="message") } })
observeEvent(input$btn_apply_outrider_filter, { df <- patient_filtered_outrider(); visible <- df[input$patient_outrider_table_rows_all, , drop=FALSE]; if(nrow(visible) > 0) { broadcast_global_filter("outrider", visible); showNotification("Global filters applied from Outrider.", type="message") } })
observeEvent(input$btn_apply_fraser_filter, { df <- patient_filtered_fraser(); visible <- df[input$patient_fraser_table_rows_all, , drop=FALSE]; if(nrow(visible) > 0) { broadcast_global_filter("fraser", visible); showNotification("Global filters applied from Fraser.", type="message") } })


# --- 3. ROW CLICK OBSERVERS ---
observeEvent(input$patient_table_output_rows_selected, {
  idx <- input$patient_table_output_rows_selected
  if(length(idx) > 0) {
    df <- isolate(selected_data_patient())
    sel_row <- df[idx[length(idx)], , drop=FALSE]
    selected_patient_variant(sel_row)
    selected_patient_parallel_variant(NULL); selected_patient_fraser_variant(NULL)
    selected_patient_ase_row(NULL); selected_patient_counts_row(NULL); selected_patient_outrider_row(NULL); selected_patient_outlier_bridge(NULL)
    
    sym_col <- grep("^(SYMBOL|hgnc_symbol|symbol|hgncSymbol)$", colnames(sel_row), value = TRUE, ignore.case = TRUE)[1]
    if(!is.na(sym_col)) selected_symbol_all(as.character(sel_row[[sym_col]])) else selected_symbol_all(NULL)
  }
}, ignoreNULL = TRUE)

observeEvent(input$patient_parallel_table_rows_selected, {
  idx <- input$patient_parallel_table_rows_selected
  if(length(idx) > 0) {
    df <- isolate(patient_filtered_parallel_data())
    sel_row <- df[idx[length(idx)], , drop=FALSE]
    selected_patient_parallel_variant(sel_row)
    selected_patient_variant(NULL); selected_patient_fraser_variant(NULL)
    selected_patient_ase_row(NULL); selected_patient_counts_row(NULL); selected_patient_outrider_row(NULL); selected_patient_outlier_bridge(NULL)
    
    sym_col <- grep("^(SYMBOL|hgnc_symbol|symbol|hgncSymbol)$", colnames(sel_row), value = TRUE, ignore.case = TRUE)[1]
    if(!is.na(sym_col)) selected_symbol_all(as.character(sel_row[[sym_col]])) else selected_symbol_all(NULL)
  }
}, ignoreNULL = TRUE)

observeEvent(input$patient_ase_table_rows_selected, {
  idx <- input$patient_ase_table_rows_selected
  if(length(idx) > 0) {
    df <- isolate(patient_filtered_ase())
    sel_row <- df[idx[length(idx)], , drop=FALSE]
    selected_patient_ase_row(sel_row); selected_patient_outlier_bridge(sel_row)
    selected_patient_variant(NULL); selected_patient_parallel_variant(NULL); selected_patient_fraser_variant(NULL)
    selected_patient_counts_row(NULL); selected_patient_outrider_row(NULL)
    
    sym_col <- grep("^(SYMBOL|hgnc_symbol|symbol|hgncSymbol)$", colnames(sel_row), value = TRUE, ignore.case = TRUE)[1]
    if(!is.na(sym_col)) selected_symbol_all(as.character(sel_row[[sym_col]])) else selected_symbol_all(NULL)
  }
}, ignoreNULL = TRUE)

observeEvent(input$patient_outrider_table_rows_selected, {
  idx <- input$patient_outrider_table_rows_selected
  if(length(idx) > 0) {
    df <- isolate(patient_filtered_outrider())
    sel_row <- df[idx[length(idx)], , drop=FALSE]
    selected_patient_outrider_row(sel_row); selected_patient_outlier_bridge(sel_row)
    selected_patient_variant(NULL); selected_patient_parallel_variant(NULL); selected_patient_fraser_variant(NULL)
    selected_patient_ase_row(NULL); selected_patient_counts_row(NULL)
    
    sym_col <- grep("^(SYMBOL|hgnc_symbol|symbol|hgncSymbol)$", colnames(sel_row), value = TRUE, ignore.case = TRUE)[1]
    if(!is.na(sym_col)) selected_symbol_all(as.character(sel_row[[sym_col]])) else selected_symbol_all(NULL)
  }
}, ignoreNULL = TRUE)

observeEvent(input$patient_fraser_table_rows_selected, {
  idx <- input$patient_fraser_table_rows_selected
  if(length(idx) > 0) {
    df <- isolate(patient_filtered_fraser())
    sel_row <- df[idx[length(idx)], , drop=FALSE]
    selected_patient_fraser_variant(sel_row)
    selected_patient_variant(NULL); selected_patient_parallel_variant(NULL)
    selected_patient_ase_row(NULL); selected_patient_counts_row(NULL); selected_patient_outrider_row(NULL); selected_patient_outlier_bridge(NULL)
    
    sym_col <- grep("^(SYMBOL|hgnc_symbol|symbol|hgncSymbol)$", colnames(sel_row), value = TRUE, ignore.case = TRUE)[1]
    if(!is.na(sym_col)) selected_symbol_all(as.character(sel_row[[sym_col]])) else selected_symbol_all(NULL)
  }
}, ignoreNULL = TRUE)

observeEvent(input$patient_counts_table_rows_selected, {
  idx <- input$patient_counts_table_rows_selected
  if(length(idx) > 0) {
    df <- isolate(patient_filtered_counts())
    sel_row <- df[idx[length(idx)], , drop=FALSE]
    selected_patient_counts_row(sel_row); selected_patient_outlier_bridge(sel_row)
    selected_patient_variant(NULL); selected_patient_parallel_variant(NULL); selected_patient_fraser_variant(NULL)
    selected_patient_ase_row(NULL); selected_patient_outrider_row(NULL)
    
    sym_col <- grep("^(SYMBOL|hgnc_symbol|symbol|hgncSymbol)$", colnames(sel_row), value = TRUE, ignore.case = TRUE)[1]
    if(!is.na(sym_col)) selected_symbol_all(as.character(sel_row[[sym_col]])) else selected_symbol_all(NULL)
  }
}, ignoreNULL = TRUE)

# --- 4. CLEAR TRIGGER ---
observeEvent(clear_trigger(), {
  req(clear_trigger())
  selected_patient_variant(NULL); selected_patient_parallel_variant(NULL); selected_patient_fraser_variant(NULL)
  selected_patient_ase_row(NULL); selected_patient_counts_row(NULL); selected_patient_outrider_row(NULL); selected_patient_outlier_bridge(NULL)
  selected_symbol_all(NULL)
  
  p_cf_main(NULL); p_cf_parallel(NULL); p_cf_ase(NULL); p_cf_counts(NULL); p_cf_outrider(NULL); p_cf_fraser(NULL)
  p_sync_main(NULL); p_sync_parallel(NULL) 
  
  last_selected_var(NULL); last_selected_fra(NULL)
  current_igv_state_var(NULL); current_igv_state_fra(NULL)
  
  dataTableProxy('patient_table_output') %>% selectRows(NULL) %>% clearSearch()
  dataTableProxy('patient_parallel_table') %>% selectRows(NULL) %>% clearSearch()
  dataTableProxy('patient_ase_table') %>% selectRows(NULL) %>% clearSearch()
  dataTableProxy('patient_counts_table') %>% selectRows(NULL) %>% clearSearch()
  dataTableProxy('patient_outrider_table') %>% selectRows(NULL) %>% clearSearch()
  dataTableProxy('patient_fraser_table') %>% selectRows(NULL) %>% clearSearch()
  
  if (length(loaded_igv_tracks_var()) > 0) session$sendCustomMessage("force-remove-igv-tracks", list(id = session$ns("igvShiny_variants"), trackNames = as.list(loaded_igv_tracks_var())))
  if (length(loaded_igv_tracks_fra()) > 0) session$sendCustomMessage("force-remove-igv-tracks", list(id = session$ns("igvShiny_fraser"), trackNames = as.list(loaded_igv_tracks_fra())))
  
  loaded_igv_tracks_var(character(0)); loaded_igv_tracks_fra(character(0))
  session$sendCustomMessage("igv-nav", list(id = session$ns("igvShiny_variants"), locus = "BRCA1"))
  session$sendCustomMessage("igv-nav", list(id = session$ns("igvShiny_fraser"), locus = "BRCA1"))
}, ignoreInit = TRUE)

# --- 5. IGV TRACK LOADING LOGIC ---
observeEvent(list(selected_patient_variant(), selected_patient_parallel_variant(), input$selected_patient), {
  req(igv_ready()) 
  var <- selected_patient_variant()
  if (is.null(var)) var <- selected_patient_parallel_variant()
  if (is.null(var)) var <- last_selected_var()
  
  current_state <- list(var = var, pids = input$selected_patient)
  if (identical(current_igv_state_var(), current_state)) return()
  current_igv_state_var(current_state) 
  
  igv_id <- session$ns("igvShiny_variants")
  if (length(loaded_igv_tracks_var()) > 0) session$sendCustomMessage("force-remove-igv-tracks", list(id = igv_id, trackNames = as.list(loaded_igv_tracks_var())))
  tryCatch({ removeUserAddedTracks(session, igv_id) }, error = function(e) {})
  loaded_igv_tracks_var(character(0))
  
  if (length(input$selected_patient) == 0) return()
  if (!is.null(var)) last_selected_var(var) 
  
  if (!is.null(var)) {
    cp <- suppressWarnings(as.numeric(var$pos)); chr <- as.character(var$chr)
    locus_str <- paste0("chr", chr, ":", cp-100, "-", cp+100)
    slice_chr <- chr; slice_start <- cp - 5000; slice_end <- cp + 5000
  } else { locus_str <- "BRCA1"; slice_chr <- "17"; slice_start <- 43044295; slice_end <- 43125483 }
  
  session$sendCustomMessage("igv-nav", list(id = igv_id, locus = locus_str))
  ts <- paste0(as.integer(Sys.time()), "_", sample(1:10000, 1)); new_tracks <- character(0)
  
  for (pid in input$selected_patient) {
    r_bam <- find_bam_file("R_RMDUP_BAM", paste0(pid, ".bam"))
    if (!is.null(r_bam)) { 
      slice_name <- paste0(session$token, "_", pid, "_R_var_", ts, ".bam")
      sl <- file.path(temp_bam_dir, slice_name)
      Rsamtools::filterBam(r_bam, sl, param=Rsamtools::ScanBamParam(which=GenomicRanges::GRanges(slice_chr, IRanges::IRanges(slice_start, slice_end))), indexDestination=TRUE)
      track_name <- paste0("RNA: ", pid, " [", substr(ts,1,4), "]")
      loadBamTrackFromURL(session, igv_id, track_name, paste0("tmp_bam_slices/", slice_name), paste0("tmp_bam_slices/", slice_name, ".bai"), TRUE, "SQUISHED") 
      new_tracks <- c(new_tracks, track_name)
    }
    if (exists("rna_to_dna") && !is.null(rna_to_dna)) {
      did <- rna_to_dna$dna_id[rna_to_dna$rna_id == pid]
      if (length(did) > 0) {
        d_bam <- find_bam_file("D_RMDUP_BAM", paste0(did[1], ".bam"))
        if (!is.null(d_bam)) {
          slice_name_d <- paste0(session$token, "_", did[1], "_D_var_", ts, ".bam")
          sl_d <- file.path(temp_bam_dir, slice_name_d)
          Rsamtools::filterBam(d_bam, sl_d, param=Rsamtools::ScanBamParam(which=GenomicRanges::GRanges(slice_chr, IRanges::IRanges(slice_start, slice_end))), indexDestination=TRUE)
          track_name_d <- paste0("DNA: ", did[1], " [", substr(ts,1,4), "]")
          loadBamTrackFromURL(session, igv_id, track_name_d, paste0("tmp_bam_slices/", slice_name_d), paste0("tmp_bam_slices/", slice_name_d, ".bai"), TRUE, "SQUISHED") 
          new_tracks <- c(new_tracks, track_name_d)
        }
      }
    }
  }
  loaded_igv_tracks_var(new_tracks)
}, ignoreInit = TRUE)

observeEvent(list(selected_patient_fraser_variant(), input$selected_patient), {
  req(igv_ready()) 
  row <- selected_patient_fraser_variant(); if (is.null(row)) row <- last_selected_fra()
  current_state <- list(var = row, pids = input$selected_patient); if (identical(current_igv_state_fra(), current_state)) return()
  current_igv_state_fra(current_state)
  
  igv_id <- session$ns("igvShiny_fraser")
  if (length(loaded_igv_tracks_fra()) > 0) session$sendCustomMessage("force-remove-igv-tracks", list(id = igv_id, trackNames = as.list(loaded_igv_tracks_fra())))
  tryCatch({ removeUserAddedTracks(session, igv_id) }, error = function(e) {})
  loaded_igv_tracks_fra(character(0))
  
  if (length(input$selected_patient) == 0) return()
  if (!is.null(row)) last_selected_fra(row) 
  
  if (!is.null(row)) {
    f_chr <- as.character(row$chr); f_start <- as.numeric(row$start); f_end <- as.numeric(row$end)
    locus_str <- paste0("chr", f_chr, ":", f_start-100, "-", f_end+100)
    slice_chr <- f_chr; slice_start <- f_start - 1000; slice_end <- f_end + 1000
  } else { locus_str <- "BRCA1"; slice_chr <- "17"; slice_start <- 43044295; slice_end <- 43125483 }
  
  session$sendCustomMessage("igv-nav", list(id = igv_id, locus = locus_str))
  ts <- paste0(as.integer(Sys.time()), "_", sample(1:10000, 1)); new_tracks <- character(0)
  
  for (pid in input$selected_patient) {
    r_bam <- find_bam_file("R_RMDUP_BAM", paste0(pid, ".bam"))
    if (!is.null(r_bam)) { 
      slice_name <- paste0(session$token, "_", pid, "_fra_", ts, ".bam")
      sl <- file.path(temp_bam_dir, slice_name)
      Rsamtools::filterBam(r_bam, sl, param=Rsamtools::ScanBamParam(which=GenomicRanges::GRanges(slice_chr, IRanges::IRanges(slice_start, slice_end))), indexDestination=TRUE)
      track_name <- paste0("RNA: ", pid, " [", substr(ts,1,4), "]")
      loadBamTrackFromURL(session, igv_id, track_name, paste0("tmp_bam_slices/", slice_name), paste0("tmp_bam_slices/", slice_name, ".bai"), TRUE, "SQUISHED") 
      new_tracks <- c(new_tracks, track_name)
    }
  }
  loaded_igv_tracks_fra(new_tracks)
}, ignoreInit = TRUE)

observeEvent(TRUE, {
  lapply(names(group_list), function(group_name) {
    group_id <- gsub(" ", "_", group_name)
    overlap <- intersect(required_columns, group_list[[group_name]])
    if (length(overlap) > 0) {
      updateCheckboxGroupInput(session, paste0("group_patient_", group_id), selected = overlap)
    }
  })
}, once = TRUE, ignoreInit = FALSE)
