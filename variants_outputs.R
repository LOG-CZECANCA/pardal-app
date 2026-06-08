# variants_outputs.R

output$table_output_highlighted <- renderDT({
  df <- variants_display_data() 
  display_df <- format_for_display(df) # APPLIES COSMETICS
  use_paging <- nrow(display_df) > 1000
  datatable(display_df, options = list(paging = use_paging, scrollX = TRUE, scrollY = if (!use_paging) "500px" else NULL, scrollCollapse = TRUE, stateSave = TRUE), filter = "top", selection = "single", rownames = FALSE)
}, server = TRUE)

output$parallel_table <- renderDT({
  df <- parallel_display_data() # USES MIRRORED DATA
  if (nrow(df) == 0) return(data.frame(Message = "No parallel data found"))
  display_df <- format_for_display(df)
  datatable(display_df, options = list(paging = nrow(display_df) > 1000, scrollX = TRUE, scrollY = if (nrow(display_df) <= 1000) "300px" else NULL, scrollCollapse = TRUE, stateSave = TRUE), filter = "top", selection = "single", class = "stripe hover row-border nowrap compact", rownames = FALSE)
})

output$ase_table <- renderDT({
  df <- filtered_ase()
  if (nrow(df) == 0) return(data.frame(Message = "No ASE data"))
  df[] <- lapply(df, as.character)
  hide_cols <- grep("^(chr_gene|start_gene|end_gene)$", colnames(df), ignore.case = TRUE, value = TRUE)
  df <- df[, setdiff(colnames(df), hide_cols), drop = FALSE]
  display_df <- format_for_display(df, sort_cols = outlier_base_columns)
  datatable(display_df, options = list(pageLength = 15, scrollX = TRUE, stateSave = TRUE), filter = "top", selection = "single", rownames = FALSE)
})

output$counts_table <- renderDT({
  df <- filtered_counts()
  if (nrow(df) == 0) return(data.frame(Message = "No Counts data"))
  hide_cols <- grep("^(chr_gene|start_gene|end_gene)$", colnames(df), ignore.case = TRUE, value = TRUE)
  df <- df[, setdiff(colnames(df), hide_cols), drop = FALSE]
  display_df <- format_for_display(df, sort_cols = outlier_base_columns)
  datatable(display_df, options = list(pageLength = 15, scrollX = TRUE, stateSave = TRUE), filter = "top", selection = "single", rownames = FALSE)
})

output$outrider_table <- renderDT({
  df <- filtered_outrider()
  if (nrow(df) == 0) return(data.frame(Message = "No Outrider data"))
  hide_cols <- grep("^(chr_gene|start_gene|end_gene)$", colnames(df), ignore.case = TRUE, value = TRUE)
  df <- df[, setdiff(colnames(df), hide_cols), drop = FALSE]
  display_df <- format_for_display(df, sort_cols = outlier_base_columns)
  datatable(display_df, options = list(pageLength = 15, scrollX = TRUE, stateSave = TRUE), filter = "top", selection = "single", rownames = FALSE)
})

output$fraser_table <- renderDT({
  df <- filtered_fraser()
  if (nrow(df) == 0) return(data.frame(Message = "No Fraser data"))
  hide_cols <- grep("^(chr_gene|start_gene|end_gene)$", colnames(df), ignore.case = TRUE, value = TRUE)
  df <- df[, setdiff(colnames(df), hide_cols), drop = FALSE]
  display_df <- format_for_display(df, sort_cols = outlier_base_columns)
  datatable(display_df, options = list(pageLength = 15, scrollX = TRUE, stateSave = TRUE), filter = "top", selection = "single", rownames = FALSE)
})

output$samples_pairs_table <- renderDT({
  file_path <- "Input/sample_pairs.txt"
  if (file.exists(file_path)) {
    df <- read.table(file_path, header = FALSE, sep = " ", stringsAsFactors = FALSE)
    if (ncol(df) >= 2) colnames(df)[1:2] <- c("RNA", "DNA") else if (ncol(df) == 1) colnames(df)[1] <- "RNA"
    datatable(df, options = list(pageLength = 50, scrollX = TRUE), rownames = FALSE, filter = "top")
  } else { data.frame(Error = "Input/sample_pairs.txt not found.") }
})

output$sample_identity_table <- renderDT({
  file_path <- "SAMPLE_IDENTITY/identical_samples.csv"
  if (file.exists(file_path)) {
    df <- read.csv(file_path, stringsAsFactors = FALSE)
    datatable(df, options = list(pageLength = 50, scrollX = TRUE), rownames = FALSE, filter = "top")
  } else { data.frame(Error = "SAMPLE_IDENTITY/identical_samples.csv not found.") }
})

output$fraser_debug_info <- renderPrint({
  sel_fraser <- selected_fraser_variant_all()
  if (is.null(sel_fraser) || nrow(sel_fraser) == 0) { cat("STATUS: No Fraser Junction Selected.\n"); return() }
  f_chr <- as.character(sel_fraser$chr[1]); f_start <- as.numeric(sel_fraser$start[1]); f_end <- as.numeric(sel_fraser$end[1])
  cat(sprintf("STATUS: Selected Junction %s:%d-%d\n", f_chr, f_start, f_end))
  if (exists("exons_gr") && !is.null(exons_gr)) {
    gtf_style <- seqlevels(exons_gr)[1]; query_chr <- f_chr
    if(grepl("^chr", gtf_style) && !grepl("^chr", f_chr)) query_chr <- paste0("chr", f_chr)
    if(!grepl("^chr", gtf_style) && grepl("^chr", f_chr)) query_chr <- gsub("^chr", "", f_chr)
    j_start <- GRanges(query_chr, IRanges(f_start, f_start)); j_end <- GRanges(query_chr, IRanges(f_end, f_end))
    cat("  > Canonical Hypothesis:\n")
    hits_s <- suppressWarnings(findOverlaps(j_start, exons_gr))
    if(length(hits_s)>0) cat(sprintf("    - Start is Exonic (Exon %d)\n", subjectHits(hits_s)[1]))
    else { left <- suppressWarnings(follow(j_start, exons_gr, ignore.strand=TRUE)); if(!is.na(left)) cat(sprintf("    - Start is Intronic (Upstream Anchor End: %d)\n", end(exons_gr[left]))) else cat("    - Start Context: Intergenic or Start of Chromosome\n") }
    hits_e <- suppressWarnings(findOverlaps(j_end, exons_gr))
    if(length(hits_e)>0) cat(sprintf("    - End is Exonic (Exon %d)\n", subjectHits(hits_e)[1]))
    else { right <- suppressWarnings(precede(j_end, exons_gr, ignore.strand=TRUE)); if(!is.na(right)) cat(sprintf("    - End is Intronic (Downstream Anchor Start: %d)\n", start(exons_gr[right]))) else cat("    - End Context: Intergenic or End of Chromosome\n") }
  } else { cat("  > (GTF not loaded)\n") }
  res <- variants_filtered_cleaned()
  cat(sprintf("  > Search Result: %d variants found matching global filters.\n", nrow(res)))
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

output$igvShiny_variants_global <- renderIgvShiny({ igvShiny(genome_options) })

output$variant_details_box_var <- renderUI({
  sel <- selected_variant_main()
  if(is.null(sel)) sel <- selected_variant_parallel()
  if(is.null(sel)) return(HTML("<span style='color: #888;'>Select a variant.</span>"))
  v_chr <- as.character(sel$chr[1]); v_pos <- as.character(sel$pos[1]); v_ref <- as.character(sel$ref[1]); v_alt <- as.character(sel$alt[1])
  match_idx <- which(as.character(variants$chr) == v_chr & as.character(variants$pos) == v_pos & as.character(variants$ref) == v_ref & as.character(variants$alt) == v_alt)
  if (length(match_idx) > 0) full_row <- variants[match_idx[1], , drop = FALSE] else full_row <- sel 
  
  rna_samples <- character(0); p_data <- filtered_parallel_data()
  if (!is.null(p_data) && nrow(p_data) > 0) {
    carrier_rows <- p_data[p_data$chr == v_chr & as.character(p_data$pos) == v_pos & p_data$ref == v_ref & p_data$alt == v_alt, ]
    rna_samples <- unique(as.character(carrier_rows$sample_id))
  }
  dna_samples <- character(0)
  if (length(rna_samples) > 0 && exists("rna_to_dna") && !is.null(rna_to_dna)) dna_samples <- unique(rna_to_dna$dna_id[rna_to_dna$rna_id %in% rna_samples])
  
  rna_sample_name <- if(length(rna_samples) > 0) paste(rna_samples, collapse = ", ") else "N/A (Filtered Out)"
  dna_sample_name <- if(length(dna_samples) > 0) paste(dna_samples, collapse = ", ") else "N/A"
  
  get_col_val <- function(patterns) {
    col <- grep(patterns, colnames(full_row), value = TRUE, ignore.case = TRUE)[1]
    if(!is.na(col)) { val <- full_row[[col]][1]; if (is.null(val) || is.na(val)) return("N/A"); return(as.character(val)) }
    return("N/A")
  }
  sym <- get_col_val("^(SYMBOL|hgnc_symbol)$"); gid <- get_col_val("^(Gene|gene_id|ENSG.*)$")
  gene_display <- if (isTRUE(gid != "N/A" && gid != "")) paste0(sym, " (", gid, ")") else sym
  refseq <- get_col_val("^(RefSeq|Refseq|refseq|NM_.*)$"); feature <- get_col_val("^(Feature|transcript|transcript_id)$")
  conseq <- get_col_val("^(Consequence|impact|effect)$"); hgvsc <- get_col_val("^(HGVSc|hgvsc)$"); clin_sig <- get_col_val("^(ClinVar_CLNSIG|CLNSIG|clinvar_clnsig|clinical_significance)$")
  
  html_content <- paste0(
    "<div style='font-family: sans-serif; font-size: 14px; line-height: 1.6;'>",
    "  <div style='margin-bottom: 10px; border-bottom: 2px solid #ccc; padding-bottom: 5px; background-color: #f1f7fc; padding: 8px; border-radius: 4px;'>",
    "    <div><strong>RNA Sample(s):</strong> <span style='color:#0055aa;'>", rna_sample_name, "</span></div>",
    "    <div><strong>DNA Sample(s):</strong> <span style='color:#0055aa;'>", dna_sample_name, "</span></div>",
    "  </div>",
    "  <div style='margin-bottom: 5px; border-bottom: 1px solid #ddd; padding-bottom: 5px;'>",
    "    <strong>Location:</strong> <span style='font-family: monospace;'>", v_chr, ":", v_pos, "</span>",
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
  var_filename <- paste0("IGV_Variant_", v_chr, "_", v_pos, ".png")
  html_content <- paste0("<div id='igv-snapshot-filename' style='display:none;'>", var_filename, "</div>", html_content)
  HTML(html_content)
})

# --- DOWNLOAD HANDLERS ---
output$dl_all_main <- downloadHandler(
  filename = function() { paste0("All_Variants_Main_", Sys.Date(), ".csv") }, 
  content = function(file) { write.csv(format_for_display(variants_display_data()), file, row.names = FALSE) }
)

output$dl_all_parallel <- downloadHandler(
  filename = function() { paste0("All_Variants_Parallel_", Sys.Date(), ".csv") }, 
  content = function(file) { 
    df <- as.data.frame(parallel_display_data())
    df[] <- lapply(df, function(x) if (is.list(x)) sapply(x, function(y) paste(as.character(y), collapse=",")) else x)
    write.csv(format_for_display(df), file, row.names = FALSE) 
  }
)

output$dl_ase <- downloadHandler(
  filename = function() { paste0("Filtered_ASE_", Sys.Date(), ".csv") },
  content = function(file) { write.csv(filtered_ase(), file, row.names = FALSE) }
)

output$dl_outrider <- downloadHandler(
  filename = function() { paste0("Filtered_Outrider_", Sys.Date(), ".csv") },
  content = function(file) { write.csv(filtered_outrider(), file, row.names = FALSE) }
)

output$dl_fraser <- downloadHandler(
  filename = function() { paste0("Filtered_Fraser_", Sys.Date(), ".csv") },
  content = function(file) { write.csv(filtered_fraser(), file, row.names = FALSE) }
)

output$dl_counts <- downloadHandler(
  filename = function() { paste0("Filtered_Counts_", Sys.Date(), ".csv") },
  content = function(file) { write.csv(filtered_counts(), file, row.names = FALSE) }
)

output$dl_all_joined <- downloadHandler(
  filename = function() { paste0("All_Variants_Joined_", Sys.Date(), ".csv") }, 
  content = function(file) { 
    main_df <- variants_display_data()
    parallel_df <- parallel_display_data()
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
    
    write_if_data(variants_display_data(), "Variants_Main.csv")
    write_if_data(parallel_display_data(), "Variants_Parallel.csv")
    write_if_data(filtered_ase(), "ASE_Outliers.csv")
    write_if_data(filtered_outrider(), "Outrider_Outliers.csv")
    write_if_data(filtered_fraser(), "Fraser_Outliers.csv")
    write_if_data(filtered_counts(), "Read_Counts.csv")
    
    if (length(files) > 0) { zip(file, files) } else { write.csv(data.frame(Message="No data generated."), file, row.names=FALSE) }
  }
)
