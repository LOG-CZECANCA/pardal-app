# variants_observers.R

if (!exists("loaded_igv_tracks_global", inherits = FALSE)) {
  loaded_igv_tracks_global <- reactiveVal(character(0))
  current_igv_state_global <- reactiveVal(NULL)
}

is_actively_searched <- function(g, c) {
  has_g <- !is.null(g) && length(g) > 0 && g != ""
  has_c <- !is.null(c) && length(c) > 0 && any(c != "")
  return(has_g || has_c)
}

observe({
  current_selections <- as.character(unlist(lapply(names(group_list), function(g) { 
    input[[paste0("group_all_", gsub(" ", "_", g))]] 
  })))
  old_order <- isolate(user_col_order_main())
  kept_cols <- old_order[old_order %in% current_selections]
  new_cols <- setdiff(current_selections, old_order)
  user_col_order_main(c(kept_cols, new_cols))
})

lapply(names(group_list), function(group_name) {
  group_id <- gsub(" ", "_", group_name)
  observeEvent(input[[paste0("toggle_all_", group_id)]], {
    new_selections <- if (isTRUE(input[[paste0("toggle_all_", group_id)]])) unname(as.character(group_list[[group_name]])) else character(0)
    updateCheckboxGroupInput(session, inputId = paste0("group_all_", group_id), selected = new_selections)
  })
})

# 🌟 GLITCH FIXED: The identical() checks prevent the tables from un-filtering each other!
observeEvent(input$table_output_highlighted_rows_all, {
  req(input$table_output_highlighted_rows_all)
  if (is_actively_searched(input$table_output_highlighted_search, input$table_output_highlighted_search_columns)) {
    df <- isolate(variants_display_data())
    idx <- input$table_output_highlighted_rows_all
    idx <- idx[!is.na(idx) & idx > 0 & idx <= nrow(df)] 
    if (length(idx) > 0) {
      new_sync <- df[idx, , drop=FALSE]
      if (!identical(new_sync, v_sync_main())) v_sync_main(new_sync)
    } else { if (!is.null(v_sync_main())) v_sync_main(NULL) }
  } else { if (!is.null(v_sync_main())) v_sync_main(NULL) }
}, ignoreNULL = FALSE, ignoreInit = TRUE)

observeEvent(input$parallel_table_rows_all, {
  req(input$parallel_table_rows_all)
  if (is_actively_searched(input$parallel_table_search, input$parallel_table_search_columns)) {
    df <- isolate(filtered_parallel_data())
    idx <- input$parallel_table_rows_all
    idx <- idx[!is.na(idx) & idx > 0 & idx <= nrow(df)] 
    if (length(idx) > 0) {
      new_sync <- df[idx, , drop=FALSE]
      if (!identical(new_sync, v_sync_parallel())) v_sync_parallel(new_sync)
    } else { if (!is.null(v_sync_parallel())) v_sync_parallel(NULL) }
  } else { if (!is.null(v_sync_parallel())) v_sync_parallel(NULL) }
}, ignoreNULL = FALSE, ignoreInit = TRUE)


broadcast_v_filter <- function(src, data) {
  payload <- list(src = src, data = data)
  
  # ONLY overwrite main/parallel if the signal came from an Outlier table!
  if(src != "main" && src != "variants") v_cf_main(payload)
  if(src != "parallel" && src != "variants") v_cf_parallel(payload)
  
  v_cf_ase(if (src == "ase") NULL else payload)
  v_cf_counts(if (src == "counts") NULL else payload)
  v_cf_outrider(if (src == "outrider") NULL else payload)
  v_cf_fraser(if (src == "fraser") NULL else payload)
}

# 🌟 SIMULTANEOUS FILTERS PRESERVED: Takes Main/Parallel without erasing your current filters!
observeEvent(input$btn_apply_global_filter, {
  df_main <- isolate(variants_display_data())
  idx_main <- input$table_output_highlighted_rows_all
  if (is.null(idx_main)) idx_main <- seq_len(nrow(df_main))
  idx_main <- idx_main[!is.na(idx_main) & idx_main > 0 & idx_main <= nrow(df_main)]
  visible_main <- df_main[idx_main, , drop=FALSE]
  
  df_par <- isolate(filtered_parallel_data())
  idx_par <- input$parallel_table_rows_all
  if (is.null(idx_par)) idx_par <- seq_len(nrow(df_par))
  idx_par <- idx_par[!is.na(idx_par) & idx_par > 0 & idx_par <= nrow(df_par)]
  visible_par <- df_par[idx_par, , drop=FALSE]
  
  if(nrow(visible_main) > 0 && nrow(visible_par) > 0) {
    # Match both tables to get the precise sample_ids for the visible rows
    matching_par <- dplyr::semi_join(visible_par, visible_main, by = intersect(colnames(visible_par), c("chr", "pos", "ref", "alt", "strand")))
    if(nrow(matching_par) > 0) {
      broadcast_v_filter("variants", matching_par)
      showNotification("Global filters safely applied. Variants tables kept intact.", type="message")
    } else { showNotification("No overlapping variants found between tables.", type="warning") }
  } else { showNotification("No variants visible to apply.", type="warning") }
})

observeEvent(input$btn_apply_ase_filter, { df <- filtered_ase(); visible <- df[input$ase_table_rows_all, , drop=FALSE]; if(nrow(visible) > 0) { broadcast_v_filter("ase", visible); showNotification("Global filters applied from ASE.", type="message") } })
observeEvent(input$btn_apply_counts_filter, { df <- filtered_counts(); visible <- df[input$counts_table_rows_all, , drop=FALSE]; if(nrow(visible) > 0) { broadcast_v_filter("counts", visible); showNotification("Global filters applied from Gene Counts.", type="message") } })
observeEvent(input$btn_apply_outrider_filter, { df <- filtered_outrider(); visible <- df[input$outrider_table_rows_all, , drop=FALSE]; if(nrow(visible) > 0) { broadcast_v_filter("outrider", visible); showNotification("Global filters applied from Outrider.", type="message") } })
observeEvent(input$btn_apply_fraser_filter, { df <- filtered_fraser(); visible <- df[input$fraser_table_rows_all, , drop=FALSE]; if(nrow(visible) > 0) { broadcast_v_filter("fraser", visible); showNotification("Global filters applied from Fraser.", type="message") } })

# ROW CLICKS (ignoreNULL = TRUE protects selections)
observeEvent(input$table_output_highlighted_rows_selected, {
  idx <- input$table_output_highlighted_rows_selected
  if(length(idx) > 0) {
    df <- isolate(variants_display_data())
    sel_row <- df[idx[length(idx)], , drop=FALSE]; selected_variant_main(sel_row)
    selected_variant_parallel(NULL); selected_fraser_variant_all(NULL); selected_ase_row(NULL); selected_outrider_row(NULL); selected_counts_row(NULL); selected_outlier_bridge(NULL)
    sym_col <- grep("^(SYMBOL|hgnc_symbol|symbol|hgncSymbol)$", colnames(sel_row), value = T, ignore.case = T)[1]; if(!is.na(sym_col)) selected_symbol_all(as.character(sel_row[[sym_col]])) else selected_symbol_all(NULL)
  }
}, ignoreNULL = TRUE)

observeEvent(input$parallel_table_rows_selected, {
  idx <- input$parallel_table_rows_selected
  if(length(idx) > 0) {
    df <- isolate(filtered_parallel_data())
    sel_row <- df[idx[length(idx)], , drop=FALSE]; selected_variant_parallel(sel_row)
    selected_variant_main(NULL); selected_fraser_variant_all(NULL); selected_ase_row(NULL); selected_outrider_row(NULL); selected_counts_row(NULL); selected_outlier_bridge(NULL)
    sym_col <- grep("^(SYMBOL|hgnc_symbol|symbol|hgncSymbol)$", colnames(sel_row), value = T, ignore.case = T)[1]; if(!is.na(sym_col)) selected_symbol_all(as.character(sel_row[[sym_col]])) else selected_symbol_all(NULL)
  }
}, ignoreNULL = TRUE)

observeEvent(input$ase_table_rows_selected, {
  idx <- input$ase_table_rows_selected
  if(length(idx) > 0) {
    df <- isolate(filtered_ase())
    sel_row <- df[idx[length(idx)], , drop=FALSE]; selected_ase_row(sel_row); selected_outlier_bridge(sel_row)
    selected_variant_main(NULL); selected_variant_parallel(NULL); selected_fraser_variant_all(NULL); selected_outrider_row(NULL); selected_counts_row(NULL)
    sym_col <- grep("^(SYMBOL|hgnc_symbol|symbol|hgncSymbol)$", colnames(sel_row), value = T, ignore.case = T)[1]; if(!is.na(sym_col)) selected_symbol_all(as.character(sel_row[[sym_col]])) else selected_symbol_all(NULL)
  }
}, ignoreNULL = TRUE)

observeEvent(input$outrider_table_rows_selected, {
  idx <- input$outrider_table_rows_selected
  if(length(idx) > 0) {
    df <- isolate(filtered_outrider())
    sel_row <- df[idx[length(idx)], , drop=FALSE]; selected_outrider_row(sel_row); selected_outlier_bridge(sel_row)
    selected_variant_main(NULL); selected_variant_parallel(NULL); selected_fraser_variant_all(NULL); selected_ase_row(NULL); selected_counts_row(NULL)
    sym_col <- grep("^(SYMBOL|hgnc_symbol|symbol|hgncSymbol)$", colnames(sel_row), value = T, ignore.case = T)[1]; if(!is.na(sym_col)) selected_symbol_all(as.character(sel_row[[sym_col]])) else selected_symbol_all(NULL)
  }
}, ignoreNULL = TRUE)

observeEvent(input$fraser_table_rows_selected, {
  idx <- input$fraser_table_rows_selected
  if(length(idx) > 0) {
    df <- isolate(filtered_fraser())
    sel_row <- df[idx[length(idx)], , drop=FALSE]; selected_fraser_variant_all(sel_row)
    selected_variant_main(NULL); selected_variant_parallel(NULL); selected_ase_row(NULL); selected_outrider_row(NULL); selected_counts_row(NULL); selected_outlier_bridge(NULL)
    sym_col <- grep("^(SYMBOL|hgnc_symbol|symbol|hgncSymbol)$", colnames(sel_row), value = T, ignore.case = T)[1]; if(!is.na(sym_col)) selected_symbol_all(as.character(sel_row[[sym_col]])) else selected_symbol_all(NULL)
  }
}, ignoreNULL = TRUE)

observeEvent(input$counts_table_rows_selected, {
  idx <- input$counts_table_rows_selected
  if(length(idx) > 0) {
    df <- isolate(filtered_counts())
    sel_row <- df[idx[length(idx)], , drop=FALSE]; selected_counts_row(sel_row); selected_outlier_bridge(sel_row)
    selected_variant_main(NULL); selected_variant_parallel(NULL); selected_fraser_variant_all(NULL); selected_ase_row(NULL); selected_outrider_row(NULL)
    sym_col <- grep("^(SYMBOL|hgnc_symbol|symbol|hgncSymbol)$", colnames(sel_row), value = T, ignore.case = T)[1]; if(!is.na(sym_col)) selected_symbol_all(as.character(sel_row[[sym_col]])) else selected_symbol_all(NULL)
  }
}, ignoreNULL = TRUE)

observeEvent(clear_trigger(), {
  req(clear_trigger())
  selected_variant_main(NULL); selected_variant_parallel(NULL); selected_fraser_variant_all(NULL)
  selected_ase_row(NULL); selected_counts_row(NULL); selected_outrider_row(NULL); selected_outlier_bridge(NULL)
  selected_gene_all(NULL); selected_sample_all(NULL); selected_symbol_all(NULL)
  
  v_cf_main(NULL); v_cf_parallel(NULL); v_cf_ase(NULL); v_cf_counts(NULL); v_cf_outrider(NULL); v_cf_fraser(NULL)
  v_sync_main(NULL); v_sync_parallel(NULL)
  
  dataTableProxy('table_output_highlighted') %>% selectRows(NULL) %>% clearSearch()
  dataTableProxy('parallel_table') %>% selectRows(NULL) %>% clearSearch()
  dataTableProxy('ase_table') %>% selectRows(NULL) %>% clearSearch()
  dataTableProxy('counts_table') %>% selectRows(NULL) %>% clearSearch()
  dataTableProxy('outrider_table') %>% selectRows(NULL) %>% clearSearch()
  dataTableProxy('fraser_table') %>% selectRows(NULL) %>% clearSearch()
  
  if (length(loaded_igv_tracks_global()) > 0) session$sendCustomMessage("force-remove-igv-tracks", list(id = session$ns("igvShiny_variants_global"), trackNames = as.list(loaded_igv_tracks_global())))
  loaded_igv_tracks_global(character(0)); current_igv_state_global(NULL)
  session$sendCustomMessage("igv-nav", list(id = session$ns("igvShiny_variants_global"), locus = "BRCA1"))
}, ignoreInit = TRUE)

observeEvent(list(selected_variant_main(), selected_variant_parallel()), {
  req(igv_ready())
  igv_id <- session$ns("igvShiny_variants_global")
  var <- selected_variant_main(); if (is.null(var)) var <- selected_variant_parallel()
  if (is.null(var)) return()
  if (identical(current_igv_state_global(), var)) return()
  current_igv_state_global(var)
  
  if (length(loaded_igv_tracks_global()) > 0) session$sendCustomMessage("force-remove-igv-tracks", list(id = igv_id, trackNames = as.list(loaded_igv_tracks_global())))
  tryCatch({ removeUserAddedTracks(session, igv_id) }, error = function(e) {})
  loaded_igv_tracks_global(character(0))
  
  cp <- suppressWarnings(as.numeric(var$pos)); chr <- as.character(var$chr)
  if (!is.na(cp)) {
    session$sendCustomMessage("igv-nav", list(id = igv_id, locus = paste0("chr", chr, ":", cp-100, "-", cp+100)))
    p_data <- isolate(filtered_parallel_data())
    matches <- p_data[p_data$chr == chr & as.character(p_data$pos) == as.character(cp) & p_data$ref == var$ref & p_data$alt == var$alt, ]
    target_patients <- unique(as.character(matches$sample_id))
    if (length(target_patients) == 0) return()
    ts <- paste0(as.integer(Sys.time()), "_", sample(1:10000, 1))
    new_tracks <- character(0)
    for (pid in target_patients) {
      r_bam <- find_bam_file("R_RMDUP_BAM", paste0(pid, ".bam"))
      if (!is.null(r_bam)) {
        sl <- file.path(temp_bam_dir, paste0(session$token, "_", pid, "_R_glob_", ts, ".bam"))
        Rsamtools::filterBam(r_bam, sl, param=Rsamtools::ScanBamParam(which=GenomicRanges::GRanges(chr, IRanges::IRanges(cp-5000, cp+5000))), indexDestination=TRUE)
        track_name <- paste0("RNA: ", pid, " [", substr(ts, 1, 4), "]")
        loadBamTrackFromURL(session, igv_id, track_name, paste0("tmp_bam_slices/", basename(sl)), paste0("tmp_bam_slices/", basename(sl), ".bai"), FALSE, "SQUISHED")
        new_tracks <- c(new_tracks, track_name)
      }
      if (exists("rna_to_dna") && !is.null(rna_to_dna)) {
        did <- rna_to_dna$dna_id[rna_to_dna$rna_id == pid]
        if (length(did) > 0) {
          d_bam <- find_bam_file("D_RMDUP_BAM", paste0(did[1], ".bam"))
          if (!is.null(d_bam)) {
            sl_d <- file.path(temp_bam_dir, paste0(session$token, "_", did[1], "_D_glob_", ts, ".bam"))
            Rsamtools::filterBam(d_bam, sl_d, param=Rsamtools::ScanBamParam(which=GenomicRanges::GRanges(chr, IRanges::IRanges(cp-5000, cp+5000))), indexDestination=TRUE)
            track_name_d <- paste0("DNA: ", did[1], " [", substr(ts, 1, 4), "]")
            loadBamTrackFromURL(session, igv_id, track_name_d, paste0("tmp_bam_slices/", basename(sl_d)), paste0("tmp_bam_slices/", basename(sl_d), ".bai"), FALSE, "SQUISHED")
            new_tracks <- c(new_tracks, track_name_d)
          }
        }
      }
    }
    loaded_igv_tracks_global(new_tracks)
  }
}, ignoreNULL = TRUE, ignoreInit = TRUE)

# 🌟 AUTO-TICKS DEFAULT COLUMNS ON STARTUP
observeEvent(TRUE, {
  lapply(names(group_list), function(group_name) {
    group_id <- gsub(" ", "_", group_name)
    overlap <- intersect(required_columns, group_list[[group_name]])
    if (length(overlap) > 0) {
      updateCheckboxGroupInput(session, paste0("group_all_", group_id), selected = overlap)
    }
  })
}, once = TRUE, ignoreInit = FALSE)