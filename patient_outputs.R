# patient_outputs.R

output$patient_table_output <- renderDT({
  df <- selected_data_patient()
  display_df <- format_for_display(df) # APPLIES COSMETICS
  use_paging <- nrow(display_df) > 1000
  datatable(display_df, options = list(paging = use_paging, scrollX = TRUE, scrollY = if (!use_paging) "500px" else NULL, scrollCollapse = TRUE, stateSave = TRUE), filter = "top", selection = "single", rownames = FALSE)
}, server = TRUE)

output$patient_parallel_table <- renderDT({
  df <- patient_parallel_display_data() # USES MIRRORED DATA
  if (nrow(df) == 0) return(data.frame(Message = "No parallel data found"))
  display_df <- format_for_display(df)
  datatable(display_df, options = list(paging = nrow(display_df) > 1000, scrollX = TRUE, scrollY = if (nrow(display_df) <= 1000) "300px" else NULL, scrollCollapse = TRUE, stateSave = TRUE), filter = "top", selection = "single", class = "stripe hover row-border nowrap compact", rownames = FALSE)
})

output$patient_ase_table <- renderDT({
  df <- patient_filtered_ase()
  if (nrow(df) == 0) return(data.frame(Message = "No ASE data"))
  df[] <- lapply(df, as.character)
  hide_cols <- grep("^(chr_gene|start_gene|end_gene)$", colnames(df), ignore.case = TRUE, value = TRUE)
  df <- df[, setdiff(colnames(df), hide_cols), drop = FALSE]
  display_df <- format_for_display(df, sort_cols = outlier_base_columns)
  datatable(display_df, options = list(pageLength = 15, scrollX = TRUE, stateSave = TRUE), filter = "top", selection = "single", rownames = FALSE)
})

output$patient_counts_table <- renderDT({
  df <- patient_filtered_counts()
  if (nrow(df) == 0) return(data.frame(Message = "No Counts data"))
  hide_cols <- grep("^(chr_gene|start_gene|end_gene)$", colnames(df), ignore.case = TRUE, value = TRUE)
  df <- df[, setdiff(colnames(df), hide_cols), drop = FALSE]
  display_df <- format_for_display(df, sort_cols = outlier_base_columns)
  datatable(display_df, options = list(pageLength = 15, scrollX = TRUE, stateSave = TRUE), filter = "top", selection = "single", rownames = FALSE)
})

output$patient_outrider_table <- renderDT({
  df <- patient_filtered_outrider()
  if (nrow(df) == 0) return(data.frame(Message = "No Outrider data for this patient"))
  hide_cols <- grep("^(chr_gene|start_gene|end_gene)$", colnames(df), ignore.case = TRUE, value = TRUE)
  df <- df[, setdiff(colnames(df), hide_cols), drop = FALSE]
  display_df <- format_for_display(df, sort_cols = outlier_base_columns)
  datatable(display_df, options = list(pageLength = 15, scrollX = TRUE, stateSave = TRUE), filter = "top", selection = "single", rownames = FALSE)
})

output$patient_fraser_table <- renderDT({
  df <- patient_filtered_fraser()
  if (nrow(df) == 0) return(data.frame(Message = "No Fraser Data found"))
  hide_cols <- grep("^(chr_gene|start_gene|end_gene)$", colnames(df), ignore.case = TRUE, value = TRUE)
  df <- df[, setdiff(colnames(df), hide_cols), drop = FALSE]
  
  sel_idx <- NULL
  saved_sel <- selected_patient_fraser_variant()
  if (!is.null(saved_sel) && nrow(saved_sel) > 0 && nrow(df) > 0) {
    matches <- which(df$chr == saved_sel$chr[1] & df$start == saved_sel$start[1] & df$end == saved_sel$end[1])
    if (length(matches) > 0) sel_idx <- matches[1]
  }
  
  display_df <- format_for_display(df, sort_cols = outlier_base_columns)
  datatable(display_df, options = list(scrollX = TRUE, stateSave = TRUE), filter = "top", selection = list(mode = "single", selected = sel_idx), rownames = FALSE)
})

output$fraser_debug_info <- renderPrint({
  sel_fraser <- selected_patient_fraser_variant()
  if (is.null(sel_fraser) || nrow(sel_fraser) == 0) { cat("STATUS: No Fraser Junction Selected.\n"); return() }
  f_chr <- as.character(sel_fraser$chr[1]); f_start <- as.numeric(sel_fraser$start[1]); f_end <- as.numeric(sel_fraser$end[1])
  cat(sprintf("STATUS: Selected Junction %s:%d-%d\n", f_chr, f_start, f_end))
  if (exists("exons_gr") && !is.null(exons_gr)) {
    gtf_style <- seqlevels(exons_gr)[1]; query_chr <- f_chr
    if(grepl("^chr", gtf_style) && !grepl("^chr", f_chr)) query_chr <- paste0("chr", f_chr)
    if(!grepl("^chr", gtf_style) && grepl("^chr", f_chr)) query_chr <- gsub("^chr", "", f_chr)
    j_start <- GRanges(query_chr, IRanges(f_start, f_start)); j_end   <- GRanges(query_chr, IRanges(f_end, f_end))
    cat("  > Canonical Hypothesis:\n")
    hits_s <- suppressWarnings(findOverlaps(j_start, exons_gr))
    if(length(hits_s)>0) { cat(sprintf("    - Start is Exonic (Exon %d)\n", subjectHits(hits_s)[1]))
    } else { left <- suppressWarnings(follow(j_start, exons_gr, ignore.strand=TRUE)); if(!is.na(left)) cat(sprintf("    - Start is Intronic (Upstream Anchor End: %d)\n", end(exons_gr[left]))) else cat("    - Start Context: Intergenic or Start of Chromosome\n") }
    hits_e <- suppressWarnings(findOverlaps(j_end, exons_gr))
    if(length(hits_e)>0) { cat(sprintf("    - End is Exonic (Exon %d)\n", subjectHits(hits_e)[1]))
    } else { right <- suppressWarnings(precede(j_end, exons_gr, ignore.strand=TRUE)); if(!is.na(right)) cat(sprintf("    - End is Intronic (Downstream Anchor Start: %d)\n", start(exons_gr[right]))) else cat("    - End Context: Intergenic or End of Chromosome\n") }
  } else { cat("  > (GTF not loaded)\n") }
  res <- patient_filtered_variants()
  cat(sprintf("  > Search Result: %d variants found in this context.\n", nrow(res)))
})

output$gtex_expression_plot <- plotly::renderPlotly({
  sym <- selected_symbol_all()
  empty_plot <- function(msg, color = "#888888") { p <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::geom_text(ggplot2::aes(x=0, y=0, label=msg), size=5, color=color); return(plotly::ggplotly(p)) }
  if (is.null(sym) || length(sym) == 0 || sym[1] == "" || sym[1] == "N/A") return(empty_plot("Select a variant or gene to view GTEx Tissue Expression"))
  gene_sym <- utils::URLencode(as.character(sym[1]))
  tryCatch({
    ref_url <- paste0("https://gtexportal.org/api/v2/reference/gene?geneId=", gene_sym)
    ref_req <- httr::GET(ref_url, httr::user_agent("Mozilla/5.0"))
    if (httr::status_code(ref_req) != 200) return(empty_plot(paste("GTEx API: Could not find gene symbol", gene_sym)))
    ref_res <- jsonlite::fromJSON(httr::content(ref_req, "text", encoding = "UTF-8"))
    if (is.null(ref_res$data) || length(ref_res$data) == 0 || nrow(ref_res$data) == 0) return(empty_plot(paste("No GTEx mapping found for gene:", gene_sym)))
    gencode_id <- ref_res$data$gencodeId[1]
    exp_url <- paste0("https://gtexportal.org/api/v2/expression/medianGeneExpression?datasetId=gtex_v8&format=json&gencodeId=", gencode_id)
    exp_req <- httr::GET(exp_url, httr::user_agent("Mozilla/5.0"))
    if (httr::status_code(exp_req) == 200) {
      exp_res <- jsonlite::fromJSON(httr::content(exp_req, "text", encoding = "UTF-8"))
      df <- exp_res$data
      if (is.null(df) || nrow(df) == 0) return(empty_plot(paste("No expression data in GTEx for:", gene_sym)))
      df <- df[order(df$median), ]; df$tissueSiteDetailId <- factor(df$tissueSiteDetailId, levels = df$tissueSiteDetailId)
      p <- ggplot2::ggplot(df, ggplot2::aes(x = tissueSiteDetailId, y = median, text = paste0("<b>Tissue:</b> ", tissueSiteDetailId, "<br><b>Median TPM:</b> ", round(median, 2)))) + ggplot2::geom_col(fill = "#FD64AB") + ggplot2::coord_flip() + ggplot2::theme_minimal() + ggplot2::labs(x = NULL, y = "Median TPM", title = paste("GTEx V8 Expression:", gene_sym)) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, color = "black", hjust = 1), axis.text.x = ggplot2::element_text(size = 10, color = "black"), axis.title = ggplot2::element_text(size = 11, face = "bold"), plot.title = ggplot2::element_text(face = "bold", size = 14), plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15))
      return(plotly::ggplotly(p, tooltip = "text"))
    } else { return(empty_plot("GTEx API Error: Could not connect to expression endpoint")) }
  }, error = function(e) { return(empty_plot(paste("Connection Error:\n", e$message), color="#cc0000")) })
})

fraser_png_data <- reactive({
  sel <- selected_patient_fraser_variant()
  if (is.null(sel) || nrow(sel) == 0) return(NULL)
  f_start   <- as.character(sel$start[1])
  f_end     <- as.character(sel$end[1])
  sample_id <- input$selected_patient[1]
  base_dir <- file.path(getwd(), "R_FRASER") 
  pattern_str <- paste0(f_start, "-", f_end, ".*", sample_id, ".*\\.pdf$")
  found_files <- list.files(base_dir, pattern = pattern_str, recursive = TRUE, full.names = TRUE)
  if (length(found_files) == 0) return(list(error = paste("Plot not found for pattern:", pattern_str)))
  pdf_path <- found_files[1]
  if (requireNamespace("pdftools", quietly = TRUE) && requireNamespace("png", quietly = TRUE)) {
    tryCatch({
      bitmap <- pdftools::pdf_render_page(pdf_path, page = 1, dpi = 150)
      tmp_png <- tempfile(fileext = ".png")
      png::writePNG(bitmap, tmp_png)
      base64_str <- base64enc::dataURI(file = tmp_png, mime = "image/png")
      unlink(tmp_png)
      return(list(src = base64_str))
    }, error = function(e) { return(list(error = paste("Conversion Error:", e$message))) })
  } else { return(list(error = "Packages 'pdftools' and 'png' are required.")) }
})

output$fraser_plot_output <- renderUI({
  data <- fraser_png_data()
  if (is.null(data)) return(HTML("<div style='color: #888; padding: 20px;'>Select a row to view the plot.</div>"))
  if (!is.null(data$error)) return(HTML(paste0("<div style='color: #8a6d3b; background: #fcf8e3; padding: 10px; border: 1px solid #faebcc;'>⚠️ ", data$error, "</div>")))
  zoom_level <- if(is.null(input$plot_zoom)) "100%" else paste0(input$plot_zoom, "%")
  tags$img(src = data$src, style = paste0("width:", zoom_level, "; max-width: none; box-shadow: 0 4px 8px rgba(0,0,0,0.1); margin: 10px;"))
})

output$igvShiny_variants <- renderIgvShiny({ igvShiny(genome_options) })
output$igvShiny_fraser <- renderIgvShiny({ igvShiny(genome_options) })

get_details_html <- function(type, session) {
  pid <- session$input$selected_patient[1]
  rna_sample_name <- if (!is.null(pid)) pid else "N/A"
  dna_sample_name <- "N/A"
  if (!is.null(pid) && exists("rna_to_dna") && !is.null(rna_to_dna)) {
    did <- rna_to_dna$dna_id[rna_to_dna$rna_id == pid]
    if (length(did) > 0) dna_sample_name <- did[1]
  }
  get_col_val <- function(row_df, patterns) {
    col <- grep(patterns, colnames(row_df), value = TRUE, ignore.case = TRUE)[1]
    if(!is.na(col)) {
      val <- row_df[[col]][1]
      if(is.list(val)) val <- paste(unlist(val), collapse = ",")
      if (is.null(val) || length(val) == 0 || is.na(val) || val == "NA") return("N/A")
      return(as.character(val))
    }
    return("N/A")
  }
  
  if (type == "fraser") {
    sel <- selected_patient_fraser_variant()
    if (is.null(sel)) return(HTML("<span style='color: #888;'>Select a splicing event to see details.</span>"))
    f_chr <- as.character(sel$chr[1]); f_start <- as.character(sel$start[1]); f_end <- as.character(sel$end[1])
    location_str <- paste0(f_chr, ":", f_start, "-", f_end)
    gene_sym <- get_col_val(sel, "^(hgncSymbol|hgnc_symbol|symbol|gene_name)$")
    impact_str <- get_col_val(sel, "^(potential_impact|impact|consequence)$")
    html_content <- paste0(
      "<div style='font-family: sans-serif; font-size: 14px; line-height: 1.6;'>",
      "  <div style='margin-bottom: 10px; border-bottom: 2px solid #ccc; padding-bottom: 5px; background-color: #f1f7fc; padding: 8px; border-radius: 4px;'>",
      "    <div><strong>RNA Sample:</strong> <span style='color:#0055aa;'>", rna_sample_name, "</span></div>",
      "    <div><strong>DNA Sample:</strong> <span style='color:#0055aa;'>", dna_sample_name, "</span></div>",
      "  </div>",
      "  <div><strong>Location:</strong> <span style='font-family: monospace;'>", location_str, "</span></div>",
      "  <div><strong>Gene:</strong> ", gene_sym, "</div>",
      "  <div><strong>Potential Impact:</strong> ", impact_str, "</div>",
      "</div>"
    )
    filename_string <- paste0("IGV_Fraser_", rna_sample_name, "_", f_chr, "_", f_start, "_", f_end, ".png")
    html_content <- paste0("<div id='igv-snapshot-filename' style='display:none;'>", filename_string, "</div>", html_content)
    return(HTML(html_content))
  } else {
    sel_main <- selected_patient_variant()
    if (is.null(sel_main)) sel_main <- selected_patient_parallel_variant()
    if (is.null(sel_main)) return(HTML("<span style='color: #888;'>Select a variant.</span>"))
    t_chr <- as.character(sel_main$chr[1]); t_pos <- as.character(sel_main$pos[1]); t_ref <- as.character(sel_main$ref[1]); t_alt <- as.character(sel_main$alt[1])
    match_idx <- which(as.character(variants$chr) == t_chr & as.character(variants$pos) == t_pos & as.character(variants$ref) == t_ref & as.character(variants$alt) == t_alt)
    if (length(match_idx) > 0) full_row <- variants[match_idx[1], , drop = FALSE] else full_row <- sel_main 
    sym <- get_col_val(full_row, "^(SYMBOL|hgnc_symbol)$"); gid <- get_col_val(full_row, "^(Gene|gene_id|ENSG.*)$")
    gene_display <- if (isTRUE(gid != "N/A" && gid != "")) paste0(sym, " (", gid, ")") else sym
    refseq <- get_col_val(full_row, "^(RefSeq|Refseq|refseq|NM_.*)$"); feature  <- get_col_val(full_row, "^(Feature|transcript|transcript_id)$")
    conseq   <- get_col_val(full_row, "^(Consequence|impact|effect)$"); hgvsc    <- get_col_val(full_row, "^(HGVSc|hgvsc)$"); clin_sig <- get_col_val(full_row, "^(ClinVar_CLNSIG|CLNSIG|clinvar_clnsig|clinical_significance)$")
    html_content <- paste0(
      "<div style='font-family: sans-serif; font-size: 14px; line-height: 1.6;'>",
      "  <div style='margin-bottom: 10px; border-bottom: 2px solid #ccc; padding-bottom: 5px; background-color: #f1f7fc; padding: 8px; border-radius: 4px;'>",
      "    <div><strong>RNA Sample:</strong> <span style='color:#0055aa;'>", rna_sample_name, "</span></div>",
      "    <div><strong>DNA Sample:</strong> <span style='color:#0055aa;'>", dna_sample_name, "</span></div>",
      "  </div>",
      "  <div style='margin-bottom: 5px; border-bottom: 1px solid #ddd; padding-bottom: 5px;'>",
      "    <strong>Location:</strong> <span style='font-family: monospace;'>", t_chr, ":", t_pos, "</span>",
      "  </div>",
      "  <table style='width: 100%; border-collapse: collapse;'>",
      "    <tr><td style='width: 110px; color: #555;'>Gene:</td><td><strong>", gene_display, "</strong></td></tr>",
      "    <tr><td style='color: #555;'>Refseq:</td><td>", refseq, "</td></tr>",
      "    <tr><td style='color: #555;'>Feature:</td><td>", feature, "</td></tr>",
      "    <tr><td style='color: #555;'>HGVSc:</td><td>", hgvsc, "</td></tr>",
      "    <tr><td style='color: #555;'>Consequence:</td><td>", conseq, "</td></tr>",
      "    <tr><td style='color: #555;'>Clinical Sig:</td><td>", clin_sig, "</td></tr>",
      "  </table>",
      "</div>"
    )
    var_filename <- paste0("IGV_Variant_", rna_sample_name, "_chr", t_chr, "_", t_pos, ".png")
    html_content <- paste0("<div id='igv-snapshot-filename' style='display:none;'>", var_filename, "</div>", html_content)
    return(HTML(html_content))
  }
}

output$variant_details_box_var <- renderUI({ get_details_html("variant", session) })
output$variant_details_box_fra <- renderUI({ get_details_html("fraser", session) })


# --- DOWNLOAD HANDLERS ---
output$dl_main <- downloadHandler(
  filename = function() { paste0("Filtered_Variants_Main_", Sys.Date(), ".csv") }, 
  content = function(file) { write.csv(format_for_display(selected_data_patient()), file, row.names = FALSE) }
)

output$dl_parallel <- downloadHandler(
  filename = function() { paste0("Filtered_Variants_Parallel_", Sys.Date(), ".csv") }, 
  content = function(file) { 
    df <- patient_parallel_display_data()
    df[] <- lapply(df, function(x) if (is.list(x)) sapply(x, function(y) paste(as.character(y), collapse=",")) else x)
    write.csv(format_for_display(df), file, row.names = FALSE) 
  }
)

output$dl_ase <- downloadHandler(
  filename = function() { paste0("Filtered_ASE_", Sys.Date(), ".csv") },
  content = function(file) { write.csv(patient_filtered_ase(), file, row.names = FALSE) }
)

output$dl_outrider <- downloadHandler(
  filename = function() { paste0("Filtered_Outrider_", Sys.Date(), ".csv") },
  content = function(file) { write.csv(patient_filtered_outrider(), file, row.names = FALSE) }
)

output$dl_fraser <- downloadHandler(
  filename = function() { paste0("Filtered_Fraser_", Sys.Date(), ".csv") },
  content = function(file) { write.csv(patient_filtered_fraser(), file, row.names = FALSE) }
)

output$dl_counts <- downloadHandler(
  filename = function() { paste0("Filtered_Counts_", Sys.Date(), ".csv") },
  content = function(file) { write.csv(patient_filtered_counts(), file, row.names = FALSE) }
)

output$dl_joined <- downloadHandler(
  filename = function() { paste0("Filtered_Variants_Joined_", Sys.Date(), ".csv") }, 
  content = function(file) { 
    main_df <- selected_data_patient()
    parallel_df <- patient_parallel_display_data()
    if (nrow(main_df) > 0 && nrow(parallel_df) > 0) { 
      parallel_df[] <- lapply(parallel_df, function(x) if (is.list(x)) sapply(x, function(y) paste(as.character(y), collapse=",")) else x)
      join_keys <- intersect(c("chr", "pos", "ref", "alt", "strand"), intersect(colnames(main_df), colnames(parallel_df)))
      joined_df <- dplyr::left_join(main_df, parallel_df, by = join_keys, suffix = c("", "_parallel")) 
    } else if (nrow(main_df) > 0) { joined_df <- main_df } else { joined_df <- data.frame(Message = "No data available to join") }
    write.csv(format_for_display(joined_df), file, row.names = FALSE) 
  }
)

# MASTER ZIP DOWNLOAD
output$dl_everything <- downloadHandler(
  filename = function() { paste0("All_Filtered_Tables_", Sys.Date(), ".zip") },
  content = function(file) {
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    files <- character()
    
    write_if_data <- function(df, fname) {
      if (!is.null(df) && nrow(df) > 0) {
        if ("data.frame" %in% class(df)) {
          df[] <- lapply(df, function(x) if(is.list(x)) sapply(x, paste, collapse=",") else x)
        }
        write.csv(df, fname, row.names = FALSE)
        files <<- c(files, fname)
      }
    }
    
    write_if_data(selected_data_patient(), "Variants_Main.csv")
    write_if_data(patient_parallel_display_data(), "Variants_Parallel.csv")
    write_if_data(patient_filtered_ase(), "ASE_Outliers.csv")
    write_if_data(patient_filtered_outrider(), "Outrider_Outliers.csv")
    write_if_data(patient_filtered_fraser(), "Fraser_Outliers.csv")
    write_if_data(patient_filtered_counts(), "Read_Counts.csv")
    
    if (length(files) > 0) { zip(file, files) } else { write.csv(data.frame(Message="No data generated."), file, row.names=FALSE) }
  }
)
