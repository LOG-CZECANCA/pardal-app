# patient_reactives.R

if (!exists("user_col_order_patient", inherits = FALSE)) {
  user_col_order_patient <- reactiveVal(character(0))
}

if (!exists("cached_gene_coords", inherits = FALSE)) {
  if (exists("exons_gr") && !is.null(exons_gr)) {
    ex_df <- as.data.frame(exons_gr)
    ex_df$seqnames_clean <- as.character(ex_df$seqnames)
    if (is.list(ex_df$gene_id) || inherits(ex_df$gene_id, "List")) {
      ex_df$gene_id_clean <- sapply(ex_df$gene_id, function(x) if(length(x)>0) as.character(x)[1] else NA)
    } else { ex_df$gene_id_clean <- as.character(ex_df$gene_id) }
    
    ex_df <- ex_df[!is.na(ex_df$gene_id_clean), ]
    cgc <- data.frame(seqnames = ex_df$seqnames_clean, gene_id = ex_df$gene_id_clean, start = as.numeric(ex_df$start), end = as.numeric(ex_df$end))
    cgc <- dplyr::summarise(dplyr::group_by(cgc, seqnames, gene_id), start_gene = min(start, na.rm = TRUE), end_gene = max(end, na.rm = TRUE), .groups = "drop")
    cgc <- as.data.frame(cgc)
    colnames(cgc)[1:2] <- c("chr_gene", "gene_id")
    cgc$chr_gene <- gsub("^chr", "", trimws(as.character(cgc$chr_gene)), ignore.case=TRUE)
    cgc$Gene_clean <- gsub("\\..*$", "", cgc$gene_id)
    cached_gene_coords <<- cgc 
  } else { cached_gene_coords <<- NULL }
}

normalize_chr <- function(x) { as.character(x); gsub("^chr", "", trimws(x), ignore.case = TRUE) }

# 🌟 PREBUILD FRASER EXTENDED BOUNDARIES
if (!exists("cached_fraser_data", inherits = FALSE)) {
  if (exists("fraser_data") && !is.null(fraser_data) && exists("exons_gr") && !is.null(exons_gr)) {
    fd <- as.data.frame(fraser_data)
    f_chr <- normalize_chr(fd$chr)
    f_st <- suppressWarnings(as.numeric(fd$start))
    f_en <- suppressWarnings(as.numeric(fd$end))
    
    bound_st <- f_st
    bound_en <- f_en
    
    gtf_style <- as.character(GenomicRanges::seqnames(exons_gr)[1])
    prefix <- if(grepl("^chr", gtf_style, ignore.case=TRUE)) "chr" else ""
    query_chr <- paste0(prefix, f_chr)
    
    st_gr <- GenomicRanges::GRanges(query_chr, IRanges::IRanges(f_st, f_st))
    en_gr <- GenomicRanges::GRanges(query_chr, IRanges::IRanges(f_en, f_en))
    
    ol_st <- suppressWarnings(GenomicRanges::findOverlaps(st_gr, exons_gr))
    if(length(ol_st) > 0) {
      st_df <- data.frame(q = S4Vectors::queryHits(ol_st), val = GenomicRanges::start(exons_gr)[S4Vectors::subjectHits(ol_st)])
      st_ag <- stats::aggregate(val ~ q, st_df, FUN = min)
      bound_st[st_ag$q] <- pmin(bound_st[st_ag$q], st_ag$val, na.rm=TRUE)
    }
    unmatched_st <- which(!seq_along(bound_st) %in% S4Vectors::queryHits(ol_st))
    if(length(unmatched_st) > 0) {
      up_st <- suppressWarnings(GenomicRanges::follow(st_gr[unmatched_st], exons_gr, ignore.strand=TRUE))
      valid <- !is.na(up_st)
      bound_st[unmatched_st[valid]] <- pmin(bound_st[unmatched_st[valid]], GenomicRanges::start(exons_gr)[up_st[valid]], na.rm=TRUE)
    }
    
    ol_en <- suppressWarnings(GenomicRanges::findOverlaps(en_gr, exons_gr))
    if(length(ol_en) > 0) {
      en_df <- data.frame(q = S4Vectors::queryHits(ol_en), val = GenomicRanges::end(exons_gr)[S4Vectors::subjectHits(ol_en)])
      en_ag <- stats::aggregate(val ~ q, en_df, FUN = max)
      bound_en[en_ag$q] <- pmax(bound_en[en_ag$q], en_ag$val, na.rm=TRUE)
    }
    unmatched_en <- which(!seq_along(bound_en) %in% S4Vectors::queryHits(ol_en))
    if(length(unmatched_en) > 0) {
      dn_en <- suppressWarnings(GenomicRanges::precede(en_gr[unmatched_en], exons_gr, ignore.strand=TRUE))
      valid <- !is.na(dn_en)
      bound_en[unmatched_en[valid]] <- pmax(bound_en[unmatched_en[valid]], GenomicRanges::end(exons_gr)[dn_en[valid]], na.rm=TRUE)
    }
    
    fd$expanded_start <- bound_st
    fd$expanded_end <- bound_en
    cached_fraser_data <<- fd
  } else if (exists("fraser_data") && !is.null(fraser_data)) {
    fd <- as.data.frame(fraser_data)
    fd$expanded_start <- fd$start
    fd$expanded_end <- fd$end
    cached_fraser_data <<- fd
  } else {
    cached_fraser_data <<- NULL
  }
}

get_col_vector_gene <- function(df) {
  c <- grep("^(Gene|gene_id|geneID)$", colnames(df), value=TRUE, ignore.case=TRUE)[1]
  if(!is.na(c)) { 
    v <- as.character(df[[c]])
    v[v=="" | v=="-" | v=="NA" | is.na(v)] <- NA
    return(gsub("\\..*$", "", v)) 
  }
  rep(NA_character_, nrow(df))
}

annotate_gene_coords <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if ("start_gene" %in% colnames(df)) return(df)
  if (!exists("cached_gene_coords") || is.null(cached_gene_coords)) return(df)
  ens_col <- grep("^(Gene|gene_id|geneID)$", colnames(df), value=TRUE, ignore.case=TRUE)[1]
  if (is.na(ens_col)) return(df)
  
  df$Gene_clean_tmp <- gsub("\\..*$", "", as.character(df[[ens_col]]))
  res <- dplyr::left_join(df, cached_gene_coords[, c("Gene_clean", "chr_gene", "start_gene", "end_gene")], by = c("Gene_clean_tmp" = "Gene_clean"))
  res$Gene_clean_tmp <- NULL
  return(res)
}

get_fraser_expanded_gr <- function(fraser_df) {
  f_chr <- normalize_chr(fraser_df$chr)
  f_st <- if("expanded_start" %in% colnames(fraser_df)) suppressWarnings(as.numeric(fraser_df$expanded_start)) else suppressWarnings(as.numeric(fraser_df$start))
  f_en <- if("expanded_end" %in% colnames(fraser_df)) suppressWarnings(as.numeric(fraser_df$expanded_end)) else suppressWarnings(as.numeric(fraser_df$end))
  f_st[is.na(f_st)] <- suppressWarnings(as.numeric(fraser_df$start))[is.na(f_st)]
  f_en[is.na(f_en)] <- suppressWarnings(as.numeric(fraser_df$end))[is.na(f_en)]
  return(GenomicRanges::GRanges(f_chr, IRanges::IRanges(f_st, f_en)))
}

filter_variants_by_fraser_gr <- function(var_df, fraser_df, p_data_ref) {
  if (is.null(fraser_df) || nrow(fraser_df) == 0 || nrow(var_df) == 0) return(var_df[0, , drop=FALSE])
  samp_col_f <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(fraser_df), value=TRUE, ignore.case=TRUE)[1]
  samp_col_p <- if(!is.null(p_data_ref)) grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(p_data_ref), value=TRUE, ignore.case=TRUE)[1] else NA
  samp_col_v <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(var_df), value=TRUE, ignore.case=TRUE)[1]
  
  f_samps_all <- c()
  if (!is.na(samp_col_f)) f_samps_all <- unique(trimws(unlist(strsplit(as.character(fraser_df[[samp_col_f]]), ","))))
  f_samps_all <- f_samps_all[f_samps_all != "" & !is.na(f_samps_all)]
  
  v_chr_arr <- normalize_chr(var_df$chr)
  f_chr_arr <- normalize_chr(fraser_df$chr)
  chr_match <- v_chr_arr %in% f_chr_arr
  
  samp_match <- rep(TRUE, nrow(var_df))
  if (length(f_samps_all) > 0) {
    if (!is.na(samp_col_v)) {
      regex_pat <- paste0("\\b(", paste(f_samps_all, collapse="|"), ")\\b")
      samp_match <- grepl(regex_pat, as.character(var_df[[samp_col_v]]))
    } else if (!is.na(samp_col_p)) {
      regex_pat <- paste0("\\b(", paste(f_samps_all, collapse="|"), ")\\b")
      p_has_samp <- grepl(regex_pat, as.character(p_data_ref[[samp_col_p]]))
      valid_p_sigs <- paste(normalize_chr(p_data_ref$chr[p_has_samp]), p_data_ref$pos[p_has_samp], p_data_ref$ref[p_has_samp], p_data_ref$alt[p_has_samp], sep="_")
      v_sig_all <- paste(v_chr_arr, var_df$pos, var_df$ref, var_df$alt, sep="_")
      samp_match <- v_sig_all %in% valid_p_sigs
    }
  }
  
  sub_var <- var_df[chr_match & samp_match, , drop=FALSE]
  if (nrow(sub_var) == 0) return(var_df[0, , drop=FALSE])
  
  target_gr <- get_fraser_expanded_gr(fraser_df)
  sub_v_chr <- normalize_chr(sub_var$chr)
  var_gr <- GenomicRanges::GRanges(sub_v_chr, IRanges::IRanges(as.numeric(sub_var$pos), as.numeric(sub_var$pos)))
  
  hits <- suppressWarnings(GenomicRanges::findOverlaps(var_gr, target_gr))
  if (length(hits) == 0) return(var_df[0, , drop=FALSE])
  
  q_hits <- S4Vectors::queryHits(hits)
  s_hits <- S4Vectors::subjectHits(hits)
  
  v_samps_list <- vector("list", nrow(sub_var))
  if (!is.na(samp_col_v)) {
    for(i in seq_len(nrow(sub_var))) v_samps_list[[i]] <- unique(trimws(unlist(strsplit(as.character(sub_var[[samp_col_v]][i]), ","))))
  } else if (!is.null(p_data_ref) && !is.na(samp_col_p)) {
    v_sig_sub <- paste(sub_v_chr, sub_var$pos, sub_var$ref, sub_var$alt, sep="_")
    p_sig_all <- paste(normalize_chr(p_data_ref$chr), p_data_ref$pos, p_data_ref$ref, p_data_ref$alt, sep="_")
    p_samps_all <- as.character(p_data_ref[[samp_col_p]])
    for(i in seq_len(nrow(sub_var))) {
      match_idx <- which(p_sig_all == v_sig_sub[i])
      v_samps_list[[i]] <- unique(trimws(unlist(strsplit(paste(p_samps_all[match_idx], collapse=","), ","))))
    }
  }
  
  f_samps_arr <- if(!is.na(samp_col_f)) as.character(fraser_df[[samp_col_f]]) else rep("", nrow(fraser_df))
  keep_final <- rep(FALSE, nrow(sub_var))
  
  for (idx in seq_along(q_hits)) {
    v_i <- q_hits[idx]
    f_i <- s_hits[idx]
    f_samp <- f_samps_arr[f_i]
    if (!is.na(f_samp) && f_samp != "") {
      f_s_vals <- unique(trimws(unlist(strsplit(f_samp, ","))))
      v_s_vals <- v_samps_list[[v_i]]
      if (any(f_s_vals %in% v_s_vals)) keep_final[v_i] <- TRUE
    } else {
      keep_final[v_i] <- TRUE
    }
  }
  return(sub_var[keep_final, , drop=FALSE])
}

filter_fraser_by_variant_gr <- function(fraser_df, var_df, p_data_ref) {
  if (is.null(var_df) || nrow(var_df) == 0 || nrow(fraser_df) == 0) return(fraser_df[0, , drop=FALSE])
  samp_col_f <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(fraser_df), value=TRUE, ignore.case=TRUE)[1]
  samp_col_p <- if(!is.null(p_data_ref)) grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(p_data_ref), value=TRUE, ignore.case=TRUE)[1] else NA
  samp_col_v <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(var_df), value=TRUE, ignore.case=TRUE)[1]
  
  v_chr_arr <- normalize_chr(var_df$chr)
  v_samps_all <- c()
  if (!is.na(samp_col_v)) {
    v_samps_all <- unique(trimws(unlist(strsplit(as.character(var_df[[samp_col_v]]), ","))))
  } else if (!is.na(samp_col_p)) {
    v_sig_all <- paste(v_chr_arr, var_df$pos, var_df$ref, var_df$alt, sep="_")
    p_sig_all <- paste(normalize_chr(p_data_ref$chr), p_data_ref$pos, p_data_ref$ref, p_data_ref$alt, sep="_")
    match_idx <- which(p_sig_all %in% v_sig_all)
    v_samps_all <- unique(trimws(unlist(strsplit(as.character(p_data_ref[[samp_col_p]][match_idx]), ","))))
  }
  v_samps_all <- v_samps_all[v_samps_all != "" & !is.na(v_samps_all)]
  
  f_chr_arr <- normalize_chr(fraser_df$chr)
  chr_match <- f_chr_arr %in% v_chr_arr
  
  samp_match <- rep(TRUE, nrow(fraser_df))
  if (length(v_samps_all) > 0 && !is.na(samp_col_f)) {
    regex_pat <- paste0("\\b(", paste(v_samps_all, collapse="|"), ")\\b")
    samp_match <- grepl(regex_pat, as.character(fraser_df[[samp_col_f]]))
  }
  
  sub_fraser <- fraser_df[chr_match & samp_match, , drop=FALSE]
  if (nrow(sub_fraser) == 0) return(fraser_df[0, , drop=FALSE])
  
  var_gr <- GenomicRanges::GRanges(v_chr_arr, IRanges::IRanges(as.numeric(var_df$pos), as.numeric(var_df$pos)))
  f_gr <- get_fraser_expanded_gr(sub_fraser)
  
  hits <- suppressWarnings(GenomicRanges::findOverlaps(f_gr, var_gr))
  if (length(hits) == 0) return(fraser_df[0, , drop=FALSE])
  
  q_hits <- S4Vectors::queryHits(hits)
  s_hits <- S4Vectors::subjectHits(hits)
  
  v_samps_list <- vector("list", nrow(var_df))
  if (!is.na(samp_col_v)) {
    for(i in seq_len(nrow(var_df))) v_samps_list[[i]] <- unique(trimws(unlist(strsplit(as.character(var_df[[samp_col_v]][i]), ","))))
  } else if (!is.null(p_data_ref) && !is.na(samp_col_p)) {
    v_sig_all <- paste(v_chr_arr, var_df$pos, var_df$ref, var_df$alt, sep="_")
    p_sig_all <- paste(normalize_chr(p_data_ref$chr), p_data_ref$pos, p_data_ref$ref, p_data_ref$alt, sep="_")
    p_samps_all <- as.character(p_data_ref[[samp_col_p]])
    for(i in seq_len(nrow(var_df))) {
      match_idx <- which(p_sig_all == v_sig_all[i])
      v_samps_list[[i]] <- unique(trimws(unlist(strsplit(paste(p_samps_all[match_idx], collapse=","), ","))))
    }
  }
  
  f_samps_sub <- if(!is.na(samp_col_f)) as.character(sub_fraser[[samp_col_f]]) else rep("", nrow(sub_fraser))
  keep_final <- rep(FALSE, nrow(sub_fraser))
  
  for (idx in seq_along(q_hits)) {
    f_i <- q_hits[idx]
    v_i <- s_hits[idx]
    f_samp <- f_samps_sub[f_i]
    if (!is.na(f_samp) && f_samp != "") {
      f_s_vals <- unique(trimws(unlist(strsplit(f_samp, ","))))
      v_s_vals <- v_samps_list[[v_i]]
      if (any(f_s_vals %in% v_s_vals)) keep_final[f_i] <- TRUE
    } else {
      keep_final[f_i] <- TRUE
    }
  }
  return(sub_fraser[keep_final, , drop=FALSE])
}

filter_outliers_direct <- function(target_df, source_df) {
  if (is.null(target_df) || nrow(target_df) == 0) return(target_df[0,,drop=FALSE])
  if (is.null(source_df) || nrow(source_df) == 0) return(target_df[0,,drop=FALSE])
  
  target_df <- annotate_gene_coords(target_df)
  source_df <- annotate_gene_coords(source_df)
  samp_col_t <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(target_df), value=TRUE, ignore.case=TRUE)[1]
  samp_col_s <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(source_df), value=TRUE, ignore.case=TRUE)[1]
  
  t_ens <- get_col_vector_gene(target_df)
  s_ens <- get_col_vector_gene(source_df)
  t_samps <- if(!is.na(samp_col_t)) as.character(target_df[[samp_col_t]]) else rep(NA, nrow(target_df))
  s_samps <- if(!is.na(samp_col_s)) as.character(source_df[[samp_col_s]]) else rep(NA, nrow(source_df))
  
  keep_target <- rep(FALSE, nrow(target_df))
  for (i in seq_len(nrow(source_df))) {
    c_s_ens <- s_ens[i]; c_s_samp <- s_samps[i]
    match_samp <- rep(TRUE, nrow(target_df))
    if (!is.na(c_s_samp) && c_s_samp != "") {
      s_s_vals <- unique(trimws(unlist(strsplit(as.character(c_s_samp), ","))))
      s_s_vals <- s_s_vals[s_s_vals != "" & !is.na(s_s_vals)]
      if (length(s_s_vals) > 0 && !all(is.na(t_samps))) {
        regex_pat <- paste0("\\b(", paste(s_s_vals, collapse="|"), ")\\b")
        match_samp <- grepl(regex_pat, t_samps)
      }
    }
    has_s_gene <- !is.na(c_s_ens)
    match_gene <- (!is.na(t_ens) & has_s_gene & (t_ens == c_s_ens))
    cur_match <- match_samp & match_gene
    cur_match[is.na(cur_match)] <- FALSE
    keep_target <- keep_target | cur_match
  }
  return(target_df[keep_target, , drop=FALSE])
}

# 🌟 HIGH PERFORMANCE VECTORIZED FUZZY FILTERS (Counts Fix + Zero Padding Fallback)
fuzzy_filter_outliers <- function(outlier_df, var_df, p_data_ref, pad = 100000, gene_only = FALSE) {
  if (is.null(outlier_df) || nrow(outlier_df) == 0) return(outlier_df[0,,drop=FALSE])
  if (is.null(var_df) || nrow(var_df) == 0) return(outlier_df[0,,drop=FALSE])
  
  outlier_df <- annotate_gene_coords(outlier_df)
  samp_col_o <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(outlier_df), value=TRUE, ignore.case=TRUE)[1]
  samp_col_v <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(var_df), value=TRUE, ignore.case=TRUE)[1]
  
  o_ens <- get_col_vector_gene(outlier_df)
  v_ens <- get_col_vector_gene(var_df)
  
  parse_samps <- function(df, col) {
    if (!is.na(col)) lapply(strsplit(as.character(df[[col]]), ","), function(x) unique(trimws(x[!is.na(x) & x != ""])))
    else replicate(nrow(df), character(0), simplify=FALSE)
  }
  
  o_samps_list <- parse_samps(outlier_df, samp_col_o)
  v_samps_list <- vector("list", nrow(var_df))
  
  if (!is.na(samp_col_v)) {
    v_samps_list <- parse_samps(var_df, samp_col_v)
  } else if (!is.null(p_data_ref)) {
    samp_col_p <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(p_data_ref), value=TRUE, ignore.case=TRUE)[1]
    if (!is.na(samp_col_p)) {
      v_sig <- paste(normalize_chr(var_df$chr), var_df$pos, var_df$ref, var_df$alt, sep="_")
      p_sig <- paste(normalize_chr(p_data_ref$chr), p_data_ref$pos, p_data_ref$ref, p_data_ref$alt, sep="_")
      sig_to_samps <- split(as.character(p_data_ref[[samp_col_p]]), p_sig)
      sig_to_samps <- lapply(sig_to_samps, function(x) unique(trimws(unlist(strsplit(x, ",")))))
      v_samps_list <- unname(sig_to_samps[v_sig])
      v_samps_list[sapply(v_samps_list, is.null)] <- list(character(0))
    }
  }
  
  keep_outlier <- rep(FALSE, nrow(outlier_df))
  
  # STEP 1: Fast Gene Match
  valid_genes <- unique(v_ens[!is.na(v_ens)])
  for (o_i in which(o_ens %in% valid_genes)) {
    if (is.na(samp_col_o)) {
      keep_outlier[o_i] <- TRUE
    } else {
      o_s_vals <- o_samps_list[[o_i]]
      if (length(o_s_vals) > 0) {
        v_idx <- which(v_ens == o_ens[o_i])
        if (any(sapply(v_samps_list[v_idx], function(x) any(x %in% o_s_vals)))) {
          keep_outlier[o_i] <- TRUE
        }
      }
    }
  }
  
  # STEP 2: Fast Coordinate Fallback
  if (!gene_only) {
    actual_pad <- 0 # 🌟 Forces true intronic boundary only
    
    v_chr <- normalize_chr(var_df$chr); v_pos <- suppressWarnings(as.numeric(var_df$pos))
    o_chr <- rep(NA_character_, nrow(outlier_df))
    if("chr_gene" %in% colnames(outlier_df)) o_chr <- normalize_chr(outlier_df$chr_gene) else if("chr" %in% colnames(outlier_df)) o_chr <- normalize_chr(outlier_df$chr)
    o_start <- rep(NA_real_, nrow(outlier_df)); o_end <- rep(NA_real_, nrow(outlier_df))
    
    if("start_gene" %in% colnames(outlier_df)) { o_start <- suppressWarnings(as.numeric(outlier_df$start_gene)) - actual_pad; o_end <- suppressWarnings(as.numeric(outlier_df$end_gene)) + actual_pad } 
    else if("start" %in% colnames(outlier_df)) { o_start <- suppressWarnings(as.numeric(outlier_df$start)) - actual_pad; o_end <- suppressWarnings(as.numeric(outlier_df$end)) + actual_pad }
    else if("pos" %in% colnames(outlier_df)) { o_start <- suppressWarnings(as.numeric(outlier_df$pos)) - actual_pad; o_end <- suppressWarnings(as.numeric(outlier_df$pos)) + actual_pad }
    
    valid_o <- !is.na(o_chr) & !is.na(o_start) & !is.na(o_end)
    valid_v <- !is.na(v_chr) & !is.na(v_pos)
    
    if (any(valid_v) && any(valid_o)) {
      v_gr <- GenomicRanges::GRanges(v_chr[valid_v], IRanges::IRanges(v_pos[valid_v], v_pos[valid_v]))
      o_gr <- GenomicRanges::GRanges(o_chr[valid_o], IRanges::IRanges(o_start[valid_o], o_end[valid_o]))
      hits <- suppressWarnings(GenomicRanges::findOverlaps(v_gr, o_gr))
      
      if (length(hits) > 0) {
        q_hits <- which(valid_v)[S4Vectors::queryHits(hits)]
        s_hits <- which(valid_o)[S4Vectors::subjectHits(hits)]
        for (idx in seq_along(q_hits)) {
          o_i <- s_hits[idx]
          if (keep_outlier[o_i]) next 
          
          v_i <- q_hits[idx]
          has_v_gene <- !is.na(v_ens[v_i])
          has_o_gene <- !is.na(o_ens[o_i])
          if (has_v_gene && has_o_gene) next 
          
          v_s_vals <- v_samps_list[[v_i]]; o_s_vals <- o_samps_list[[o_i]]
          samp_match_ok <- is.na(samp_col_o) || (length(v_s_vals) > 0 && length(o_s_vals) > 0 && any(v_s_vals %in% o_s_vals))
          if (samp_match_ok) keep_outlier[o_i] <- TRUE
        }
      }
    }
  }
  return(outlier_df[keep_outlier, , drop=FALSE])
}

fuzzy_filter_variants <- function(var_df, outlier_df, p_data_ref, pad = 100000, gene_only = FALSE) {
  if (is.null(var_df) || nrow(var_df) == 0) return(var_df[0,,drop=FALSE])
  if (is.null(outlier_df) || nrow(outlier_df) == 0) return(var_df[0,,drop=FALSE])
  
  outlier_df <- annotate_gene_coords(outlier_df)
  samp_col_o <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(outlier_df), value=TRUE, ignore.case=TRUE)[1]
  samp_col_v <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(var_df), value=TRUE, ignore.case=TRUE)[1]
  
  o_ens <- get_col_vector_gene(outlier_df)
  v_ens <- get_col_vector_gene(var_df)
  
  parse_samps <- function(df, col) {
    if (!is.na(col)) lapply(strsplit(as.character(df[[col]]), ","), function(x) unique(trimws(x[!is.na(x) & x != ""])))
    else replicate(nrow(df), character(0), simplify=FALSE)
  }
  
  o_samps_list <- parse_samps(outlier_df, samp_col_o)
  v_samps_list <- vector("list", nrow(var_df))
  
  if (!is.na(samp_col_v)) {
    v_samps_list <- parse_samps(var_df, samp_col_v)
  } else if (!is.null(p_data_ref)) {
    samp_col_p <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(p_data_ref), value=TRUE, ignore.case=TRUE)[1]
    if (!is.na(samp_col_p)) {
      v_sig <- paste(normalize_chr(var_df$chr), var_df$pos, var_df$ref, var_df$alt, sep="_")
      p_sig <- paste(normalize_chr(p_data_ref$chr), p_data_ref$pos, p_data_ref$ref, p_data_ref$alt, sep="_")
      sig_to_samps <- split(as.character(p_data_ref[[samp_col_p]]), p_sig)
      sig_to_samps <- lapply(sig_to_samps, function(x) unique(trimws(unlist(strsplit(x, ",")))))
      v_samps_list <- unname(sig_to_samps[v_sig])
      v_samps_list[sapply(v_samps_list, is.null)] <- list(character(0))
    }
  }
  
  keep_var <- rep(FALSE, nrow(var_df))
  
  valid_genes <- unique(o_ens[!is.na(o_ens)])
  for (v_i in which(v_ens %in% valid_genes)) {
    if (is.na(samp_col_o)) {
      keep_var[v_i] <- TRUE
    } else {
      v_s_vals <- v_samps_list[[v_i]]
      if (length(v_s_vals) > 0) {
        o_idx <- which(o_ens == v_ens[v_i])
        if (any(sapply(o_samps_list[o_idx], function(x) any(x %in% v_s_vals)))) {
          keep_var[v_i] <- TRUE
        }
      }
    }
  }
  
  if (!gene_only) {
    actual_pad <- 0 # 🌟 Forces true intronic boundary only
    
    v_chr <- normalize_chr(var_df$chr); v_pos <- suppressWarnings(as.numeric(var_df$pos))
    o_chr <- rep(NA_character_, nrow(outlier_df))
    if("chr_gene" %in% colnames(outlier_df)) o_chr <- normalize_chr(outlier_df$chr_gene) else if("chr" %in% colnames(outlier_df)) o_chr <- normalize_chr(outlier_df$chr)
    o_start <- rep(NA_real_, nrow(outlier_df)); o_end <- rep(NA_real_, nrow(outlier_df))
    
    if("start_gene" %in% colnames(outlier_df)) { o_start <- suppressWarnings(as.numeric(outlier_df$start_gene)) - actual_pad; o_end <- suppressWarnings(as.numeric(outlier_df$end_gene)) + actual_pad } 
    else if("start" %in% colnames(outlier_df)) { o_start <- suppressWarnings(as.numeric(outlier_df$start)) - actual_pad; o_end <- suppressWarnings(as.numeric(outlier_df$end)) + actual_pad }
    else if("pos" %in% colnames(outlier_df)) { o_start <- suppressWarnings(as.numeric(outlier_df$pos)) - actual_pad; o_end <- suppressWarnings(as.numeric(outlier_df$pos)) + actual_pad }
    
    valid_o <- !is.na(o_chr) & !is.na(o_start) & !is.na(o_end)
    valid_v <- !is.na(v_chr) & !is.na(v_pos)
    
    if (any(valid_v) && any(valid_o)) {
      v_gr <- GenomicRanges::GRanges(v_chr[valid_v], IRanges::IRanges(v_pos[valid_v], v_pos[valid_v]))
      o_gr <- GenomicRanges::GRanges(o_chr[valid_o], IRanges::IRanges(o_start[valid_o], o_end[valid_o]))
      hits <- suppressWarnings(GenomicRanges::findOverlaps(v_gr, o_gr))
      
      if (length(hits) > 0) {
        q_hits <- which(valid_v)[S4Vectors::queryHits(hits)]
        s_hits <- which(valid_o)[S4Vectors::subjectHits(hits)]
        for (idx in seq_along(q_hits)) {
          v_i <- q_hits[idx]
          if (keep_var[v_i]) next
          
          o_i <- s_hits[idx]
          has_v_gene <- !is.na(v_ens[v_i])
          has_o_gene <- !is.na(o_ens[o_i])
          if (has_v_gene && has_o_gene) next 
          
          v_s_vals <- v_samps_list[[v_i]]; o_s_vals <- o_samps_list[[o_i]]
          samp_match_ok <- is.na(samp_col_o) || (length(v_s_vals) > 0 && length(o_s_vals) > 0 && any(v_s_vals %in% o_s_vals))
          if (samp_match_ok) keep_var[v_i] <- TRUE
        }
      }
    }
  }
  return(var_df[keep_var, , drop=FALSE])
}

# --- 1. BASE DATASETS ---
# --- 1. BASE DATASETS ---
base_main_patient <- reactive({
  req(input$selected_patient)
  df <- as.data.frame(variants)
  if (isTRUE(input$filter_af) && "AF" %in% colnames(df)) { af_num <- suppressWarnings(as.numeric(as.character(df$AF))); df <- df[which(is.na(af_num) | af_num < 0.01), , drop = FALSE] }
  
  # 🌟 NEW CLINICAL SIGNIFICANCE FILTER
  if (isTRUE(input$filter_clinsig)) { 
    clin_col <- grep("^(ClinVar_CLNSIG|CLNSIG|CLIN_SIG|clinical_significance)$", colnames(df), value = TRUE, ignore.case = TRUE)[1]
    if (!is.na(clin_col)) {
      c_val <- as.character(df[[clin_col]])
      keep_mask <- !grepl("benign", c_val, ignore.case = TRUE) | 
        grepl("uncertain", c_val, ignore.case = TRUE) | 
        grepl("pathogenic", c_val, ignore.case = TRUE) | 
        is.na(c_val)
      df <- df[which(keep_mask), , drop = FALSE]
    }
  }
  
  if (isTRUE(input$filter_feature) && "Feature" %in% colnames(df)) { df <- df[which(!is.na(df$Feature) & trimws(as.character(df$Feature)) != "" & trimws(as.character(df$Feature)) != "NA"), , drop = FALSE] }
  p_data <- parallel_data[parallel_data$sample_id %in% input$selected_patient, , drop = FALSE]
  if (isTRUE(input$filter_gq) && nrow(p_data) > 0) { d_gq <- suppressWarnings(as.numeric(as.character(p_data$dna_GQ))); r_gq <- suppressWarnings(as.numeric(as.character(p_data$rna_GQ))); p_data <- p_data[which(d_gq == 99 | r_gq == 99), , drop = FALSE] }
  if (isTRUE(input$filter_dp) && nrow(p_data) > 0) { d_dp <- suppressWarnings(as.numeric(as.character(p_data$dna_DP))); r_dp <- suppressWarnings(as.numeric(as.character(p_data$rna_DP))); p_data <- p_data[which(d_dp > 10 | r_dp > 10), , drop = FALSE] }
  if (nrow(p_data) > 0) df <- dplyr::semi_join(df, unique(p_data[, c("chr", "pos", "ref", "alt", "strand")]), by = c("chr", "pos", "ref", "alt", "strand")) else df <- df[0, , drop = FALSE]
  return(df)
})

base_parallel_patient <- reactive({
  req(input$selected_patient)
  df <- parallel_data[parallel_data$sample_id %in% input$selected_patient, , drop = FALSE]
  if (isTRUE(input$filter_gq) && nrow(df) > 0) { d_gq <- suppressWarnings(as.numeric(as.character(df$dna_GQ))); r_gq <- suppressWarnings(as.numeric(as.character(df$rna_GQ))); df <- df[which(d_gq == 99 | r_gq == 99), , drop = FALSE] }
  if (isTRUE(input$filter_dp) && nrow(df) > 0) { d_dp <- suppressWarnings(as.numeric(as.character(df$dna_DP))); r_dp <- suppressWarnings(as.numeric(as.character(df$rna_DP))); df <- df[which(d_dp > 10 | r_dp > 10), , drop = FALSE] }
  if (isTRUE(input$filter_af) && "AF" %in% colnames(df)) { af_num <- suppressWarnings(as.numeric(as.character(df$AF))); df <- df[which(is.na(af_num) | af_num < 0.01), , drop = FALSE] }
  
  # 🌟 NEW CLINICAL SIGNIFICANCE FILTER
  if (isTRUE(input$filter_clinsig)) { 
    clin_col <- grep("^(ClinVar_CLNSIG|CLNSIG|CLIN_SIG|clinical_significance)$", colnames(df), value = TRUE, ignore.case = TRUE)[1]
    if (!is.na(clin_col)) {
      c_val <- as.character(df[[clin_col]])
      keep_mask <- !grepl("benign", c_val, ignore.case = TRUE) | 
        grepl("uncertain", c_val, ignore.case = TRUE) | 
        grepl("pathogenic", c_val, ignore.case = TRUE) | 
        is.na(c_val)
      df <- df[which(keep_mask), , drop = FALSE]
    }
  }
  
  if (isTRUE(input$filter_feature) && "Feature" %in% colnames(df)) { df <- df[which(!is.na(df$Feature) & trimws(as.character(df$Feature)) != "" & trimws(as.character(df$Feature)) != "NA"), , drop = FALSE] }
  df <- dplyr::semi_join(df, base_main_patient(), by = c("chr", "pos", "ref", "alt", "strand"))
  return(df)
})

base_ase_patient <- reactive({
  req(input$selected_patient)
  if (!exists("ase_data") || is.null(ase_data)) return(data.frame(Message = "ASE data not loaded"))
  df <- ase_data
  samp_col <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(df), value = TRUE, ignore.case = TRUE)[1]
  if (!is.na(samp_col)) df <- df[as.character(df[[samp_col]]) %in% input$selected_patient, , drop = FALSE]
  return(df)
})

base_counts_patient <- reactive({
  req(input$selected_patient)
  if (!exists("batch_counts") || is.null(batch_counts)) return(data.frame(Message = "Counts data not loaded"))
  df <- batch_counts
  meta_regex <- "^(gene|gene_id|geneid|hgnc_symbol|symbol|chr|start|end|strand|length)$"
  meta_cols <- grep(meta_regex, colnames(df), value = TRUE, ignore.case = TRUE)
  all_cols <- colnames(df)
  pat_cols <- c()
  for (p in input$selected_patient) {
    exact <- all_cols[tolower(all_cols) == tolower(p)]
    if (length(exact) > 0) { pat_cols <- c(pat_cols, exact) } else { pat_cols <- c(pat_cols, all_cols[grepl(p, all_cols, ignore.case = TRUE)]) }
  }
  pat_cols <- unique(pat_cols)
  final_cols <- unique(c(meta_cols, pat_cols))
  if (length(final_cols) > 0) df <- df[, final_cols, drop = FALSE]
  
  if (isTRUE(input$filter_low_counts) && length(pat_cols) > 0) {
    counts_mat <- df[, pat_cols, drop = FALSE]
    counts_mat[] <- lapply(counts_mat, function(x) suppressWarnings(as.numeric(as.character(x))))
    df <- df[which(rowSums(counts_mat > 3, na.rm = TRUE) > 0), , drop = FALSE]
  }
  return(df)
})

base_outrider_patient <- reactive({
  req(input$selected_patient)
  if (!exists("outrider_data") || is.null(outrider_data)) return(data.frame(Message = "Outrider data not loaded"))
  df <- outrider_data
  samp_col <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(df), value = TRUE, ignore.case = TRUE)[1]
  if (!is.na(samp_col)) df <- df[as.character(df[[samp_col]]) %in% input$selected_patient, , drop = FALSE]
  return(df)
})

base_fraser_patient <- reactive({
  req(input$selected_patient)
  if (!exists("cached_fraser_data") || is.null(cached_fraser_data)) return(data.frame(Message = "Fraser data not loaded"))
  df <- cached_fraser_data
  samp_col <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(df), value = TRUE, ignore.case = TRUE)[1]
  if (!is.na(samp_col)) df <- df[as.character(df[[samp_col]]) %in% input$selected_patient, , drop = FALSE]
  return(df)
})

# --- 2. FILTERED TABLES ---
patient_filtered_variants <- reactive({
  df <- base_main_patient()
  if (!is.null(selected_patient_variant())) return(dplyr::semi_join(df, as.data.frame(selected_patient_variant()), by=c("chr", "pos", "ref", "alt", "strand")))
  
  gb <- p_cf_main() 
  if (!is.null(gb) && gb$src == "main") return(df)
  
  p_df <- selected_patient_parallel_variant(); if (is.null(p_df) && !is.null(gb) && gb$src == "parallel") p_df <- gb$data
  if (!is.null(p_df) && nrow(p_df) > 0) return(dplyr::semi_join(df, as.data.frame(p_df), by=c("chr", "pos", "ref", "alt", "strand")))
  
  f_df <- selected_patient_fraser_variant(); if (is.null(f_df) && !is.null(gb) && gb$src == "fraser") f_df <- gb$data
  if (!is.null(f_df) && nrow(f_df) > 0) return(filter_variants_by_fraser_gr(df, f_df, base_parallel_patient()))
  
  o_df <- selected_patient_ase_row(); if (is.null(o_df)) o_df <- selected_patient_outrider_row(); if (is.null(o_df)) o_df <- selected_patient_counts_row()
  if (is.null(o_df) && !is.null(gb) && gb$src %in% c("ase", "counts", "outrider")) o_df <- gb$data
  if (!is.null(o_df) && nrow(o_df) > 0) {
    if(!is.null(gb) && gb$src == "counts") return(fuzzy_filter_variants(df, o_df, base_parallel_patient(), pad=100000, gene_only=TRUE))
    return(fuzzy_filter_variants(df, o_df, base_parallel_patient(), pad=100000))
  }
  
  sm <- p_sync_parallel()
  if (!is.null(sm) && nrow(sm) > 0) return(dplyr::semi_join(df, sm, by=c("chr", "pos", "ref", "alt", "strand")))
  return(df)
})

patient_filtered_parallel_data <- reactive({
  df <- base_parallel_patient()
  if (!is.null(selected_patient_parallel_variant())) return(dplyr::semi_join(df, as.data.frame(selected_patient_parallel_variant()), by=c("chr", "pos", "ref", "alt", "strand")))
  
  gb <- p_cf_parallel(); if (!is.null(gb) && gb$src == "parallel") return(df)
  
  m_df <- selected_patient_variant(); if (is.null(m_df) && !is.null(gb) && gb$src == "main") m_df <- gb$data
  if (!is.null(m_df) && nrow(m_df) > 0) return(dplyr::semi_join(df, as.data.frame(m_df), by=c("chr", "pos", "ref", "alt", "strand")))
  
  f_df <- selected_patient_fraser_variant(); if (is.null(f_df) && !is.null(gb) && gb$src == "fraser") f_df <- gb$data
  if (!is.null(f_df) && nrow(f_df) > 0) {
    v_sub <- filter_variants_by_fraser_gr(base_main_patient(), f_df, df)
    return(dplyr::semi_join(df, v_sub, by=c("chr","pos","ref","alt","strand")))
  }
  
  o_df <- NULL
  if (!is.null(selected_patient_ase_row())) o_df <- selected_patient_ase_row()
  else if (!is.null(selected_patient_outrider_row())) o_df <- selected_patient_outrider_row()
  else if (!is.null(selected_patient_counts_row())) o_df <- selected_patient_counts_row()
  else if (!is.null(selected_patient_outlier_bridge())) o_df <- selected_patient_outlier_bridge()
  
  if (is.null(o_df) && !is.null(gb) && gb$src %in% c("ase", "counts", "outrider")) o_df <- gb$data
  if (!is.null(o_df) && nrow(o_df) > 0) {
    if(!is.null(gb) && gb$src == "counts") v_sub <- fuzzy_filter_variants(base_main_patient(), o_df, df, pad=100000, gene_only=TRUE)
    else v_sub <- fuzzy_filter_variants(base_main_patient(), o_df, df, pad=100000)
    return(dplyr::semi_join(df, v_sub, by=c("chr","pos","ref","alt","strand")))
  }
  
  sm <- p_sync_main()
  if (!is.null(sm) && nrow(sm) > 0) return(dplyr::semi_join(df, sm, by=c("chr", "pos", "ref", "alt", "strand")))
  return(df)
})

# --- 4. DOWNSTREAM OUTLIER TABLES ---
patient_filtered_ase <- reactive({
  if (!exists("ase_data") || is.null(ase_data)) return(NULL)
  df <- annotate_gene_coords(base_ase_patient())
  if (!is.null(selected_patient_ase_row())) return(dplyr::semi_join(df, as.data.frame(selected_patient_ase_row()), by=intersect(colnames(df), colnames(selected_patient_ase_row()))))
  gb <- p_cf_ase(); if (!is.null(gb) && gb$src == "ase") return(df)
  
  var_df <- selected_patient_variant(); if (is.null(var_df)) var_df <- selected_patient_parallel_variant()
  if (is.null(var_df) && !is.null(gb) && gb$src %in% c("main", "parallel", "variants")) var_df <- gb$data
  if (!is.null(var_df) && nrow(var_df) > 0) return(fuzzy_filter_outliers(df, var_df, base_parallel_patient(), pad=100000))
  
  out_df <- selected_patient_outrider_row(); if (is.null(out_df)) out_df <- selected_patient_counts_row(); if (is.null(out_df)) out_df <- selected_patient_fraser_variant()
  if (is.null(out_df) && !is.null(selected_patient_outlier_bridge())) out_df <- selected_patient_outlier_bridge()
  if (is.null(out_df) && !is.null(gb) && gb$src %in% c("counts", "outrider", "fraser")) out_df <- gb$data
  if (!is.null(out_df) && nrow(out_df) > 0) {
    if(!is.null(gb) && gb$src == "fraser") {
      v_bridge <- filter_variants_by_fraser_gr(base_main_patient(), out_df, base_parallel_patient())
      return(fuzzy_filter_outliers(df, v_bridge, base_parallel_patient(), pad=100000))
    }
    return(filter_outliers_direct(df, out_df))
  }
  return(df)
})

patient_filtered_outrider <- reactive({
  if (!exists("outrider_data") || is.null(outrider_data)) return(NULL)
  df <- annotate_gene_coords(base_outrider_patient())
  if (!is.null(selected_patient_outrider_row())) return(dplyr::semi_join(df, as.data.frame(selected_patient_outrider_row()), by=intersect(colnames(df), colnames(selected_patient_outrider_row()))))
  gb <- p_cf_outrider(); if (!is.null(gb) && gb$src == "outrider") return(df)
  
  var_df <- selected_patient_variant(); if (is.null(var_df)) var_df <- selected_patient_parallel_variant()
  if (is.null(var_df) && !is.null(gb) && gb$src %in% c("main", "parallel", "variants")) var_df <- gb$data
  if (!is.null(var_df) && nrow(var_df) > 0) return(fuzzy_filter_outliers(df, var_df, base_parallel_patient(), pad=100000))
  
  out_df <- selected_patient_ase_row(); if (is.null(out_df)) out_df <- selected_patient_counts_row(); if (is.null(out_df)) out_df <- selected_patient_fraser_variant()
  if (is.null(out_df) && !is.null(selected_patient_outlier_bridge())) out_df <- selected_patient_outlier_bridge()
  if (is.null(out_df) && !is.null(gb) && gb$src %in% c("counts", "ase", "fraser")) out_df <- gb$data
  if (!is.null(out_df) && nrow(out_df) > 0) {
    if(!is.null(gb) && gb$src == "fraser") {
      v_bridge <- filter_variants_by_fraser_gr(base_main_patient(), out_df, base_parallel_patient())
      return(fuzzy_filter_outliers(df, v_bridge, base_parallel_patient(), pad=100000))
    }
    return(filter_outliers_direct(df, out_df))
  }
  return(df)
})

patient_filtered_fraser <- reactive({
  if (!exists("fraser_data") || is.null(fraser_data)) return(NULL)
  df <- annotate_gene_coords(base_fraser_patient())
  if (!is.null(selected_patient_fraser_variant())) return(dplyr::semi_join(df, as.data.frame(selected_patient_fraser_variant()), by=intersect(colnames(df), colnames(selected_patient_fraser_variant()))))
  gb <- p_cf_fraser(); if (!is.null(gb) && gb$src == "fraser") return(df)
  
  var_df <- selected_patient_parallel_variant(); if (is.null(var_df)) var_df <- selected_patient_variant()
  if (is.null(var_df) && !is.null(gb) && gb$src %in% c("main", "parallel", "variants")) var_df <- gb$data
  if (!is.null(var_df) && nrow(var_df) > 0) return(filter_fraser_by_variant_gr(df, var_df, base_parallel_patient())) 
  
  out_df <- selected_patient_outrider_row(); if (is.null(out_df)) out_df <- selected_patient_counts_row(); if (is.null(out_df)) out_df <- selected_patient_ase_row()
  if (is.null(out_df) && !is.null(selected_patient_outlier_bridge())) out_df <- selected_patient_outlier_bridge()
  if (is.null(out_df) && !is.null(gb) && gb$src %in% c("counts", "outrider", "ase")) out_df <- gb$data
  if (!is.null(out_df) && nrow(out_df) > 0) {
    if(!is.null(gb) && gb$src == "counts") v_bridge <- fuzzy_filter_variants(base_main_patient(), out_df, base_parallel_patient(), pad=100000, gene_only=TRUE)
    else v_bridge <- fuzzy_filter_variants(base_main_patient(), out_df, base_parallel_patient(), pad=100000)
    return(filter_fraser_by_variant_gr(df, v_bridge, base_parallel_patient()))
  }
  return(df)
})

# 🌟 HIGH SPEED COUNTS REACTIVE
# 🌟 HIGH SPEED COUNTS REACTIVE (Coordinate Fallback Enabled!)
patient_filtered_counts <- reactive({
  if (!exists("batch_counts") || is.null(batch_counts)) return(NULL)
  df <- annotate_gene_coords(base_counts_patient())
  if (!is.null(selected_patient_counts_row())) return(dplyr::semi_join(df, as.data.frame(selected_patient_counts_row()), by=intersect(colnames(df), colnames(selected_patient_counts_row()))))
  gb <- p_cf_counts(); if (!is.null(gb) && gb$src == "counts") return(df)
  
  var_df <- selected_patient_variant(); if (is.null(var_df)) var_df <- selected_patient_parallel_variant()
  if (is.null(var_df) && !is.null(gb) && gb$src %in% c("main", "parallel", "variants")) var_df <- gb$data
  
  if (!is.null(var_df) && nrow(var_df) > 0) {
    # REMOVED gene_only=TRUE
    df <- fuzzy_filter_outliers(df, var_df, base_parallel_patient(), pad=100000, gene_only=FALSE)
    
    p_data <- base_parallel_patient()
    samp_col <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(p_data), value = TRUE, ignore.case = TRUE)[1]
    
    samps <- c()
    if (!is.na(samp_col)) {
      p_sigs <- paste(normalize_chr(p_data$chr), p_data$pos, p_data$ref, p_data$alt, sep="_")
      v_sigs <- paste(normalize_chr(var_df$chr), var_df$pos, var_df$ref, var_df$alt, sep="_")
      match_idx <- which(p_sigs %in% v_sigs)
      samps <- as.character(p_data[[samp_col]][match_idx])
    }
    
    if (length(samps)>0) {
      s_split <- unique(trimws(unlist(strsplit(as.character(samps), ","))))
      s_split <- s_split[!is.na(s_split) & s_split != ""]
      meta_cols <- grep("^(gene|gene_id|geneid|hgnc_symbol|symbol|chr|start|end|strand|length|chr_gene|start_gene|end_gene)$", colnames(df), value=T, ignore.case=T)
      pat_cols <- unlist(lapply(s_split, function(p) { colnames(df)[tolower(colnames(df)) == tolower(p)] }))
      if (length(c(meta_cols, pat_cols))>0) df <- df[, unique(c(meta_cols, pat_cols)), drop=FALSE]
    }
    return(df)
  }
  
  out_df <- selected_patient_ase_row(); if (is.null(out_df)) out_df <- selected_patient_outrider_row(); if (is.null(out_df)) out_df <- selected_patient_fraser_variant()
  if (is.null(out_df) && !is.null(selected_patient_outlier_bridge())) out_df <- selected_patient_outlier_bridge()
  if (is.null(out_df) && !is.null(gb) && gb$src %in% c("ase", "outrider", "fraser")) out_df <- gb$data
  if (!is.null(out_df) && nrow(out_df) > 0) {
    if(!is.null(gb) && gb$src == "fraser") {
      v_bridge <- filter_variants_by_fraser_gr(base_main_patient(), out_df, base_parallel_patient())
      df <- fuzzy_filter_outliers(df, v_bridge, base_parallel_patient(), pad=100000, gene_only=FALSE)
    } else {
      df <- filter_outliers_direct(df, out_df)
    }
    
    samp_col_o <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(out_df), value=TRUE, ignore.case=TRUE)[1]
    if (!is.na(samp_col_o)) {
      s_split <- unique(trimws(unlist(strsplit(as.character(out_df[[samp_col_o]]), ","))))
      s_split <- s_split[!is.na(s_split) & s_split != ""]
      if (length(s_split) > 0) {
        meta_cols <- grep("^(gene|gene_id|geneid|hgnc_symbol|symbol|chr|start|end|strand|length|chr_gene|start_gene|end_gene)$", colnames(df), value=T, ignore.case=T)
        pat_cols <- unlist(lapply(s_split, function(p) { colnames(df)[tolower(colnames(df)) == tolower(p)] }))
        if (length(c(meta_cols, pat_cols))>0) df <- df[, unique(c(meta_cols, pat_cols)), drop=FALSE]
      }
    }
  }
  return(df)
})

# --- 5. DISPLAY ASSEMBLY ---
selected_data_patient <- reactive({
  anchor_cols <- main_anchor_columns 
  user_selected_cols <- user_col_order_patient()
  
  if (length(user_selected_cols) == 0) {
    dynamic_cols <- setdiff(required_columns, anchor_cols)
  } else {
    dynamic_cols <- setdiff(user_selected_cols, anchor_cols)
  }
  final_order <- unique(c(anchor_cols, dynamic_cols))
  
  df <- patient_filtered_variants()
  if (is.null(df) || nrow(df) == 0) return(data.frame(Message = "No variants matching filters"))
  
  if ("Genomic_Context" %in% final_order && exists("exons_gr") && !is.null(exons_gr) && exists("introns_gr")) {
    v_chr <- normalize_chr(df$chr); v_pos <- as.numeric(df$pos)
    gtf_chr <- as.character(GenomicRanges::seqnames(exons_gr)[1])
    var_gr <- GenomicRanges::GRanges(paste0(if(grepl("^chr", gtf_chr)) "chr" else "", v_chr), IRanges::IRanges(v_pos, v_pos))
    h_exon <- suppressWarnings(GenomicRanges::countOverlaps(var_gr, exons_gr) > 0)
    h_intron <- suppressWarnings(GenomicRanges::countOverlaps(var_gr, introns_gr) > 0)
    df$Genomic_Context <- ifelse(h_exon & h_intron, "Exon/Intron Boundary", 
                                 ifelse(h_exon, "Exonic", 
                                        ifelse(h_intron, "Intronic", "Intergenic")))
  }
  
  df <- df[, intersect(final_order, colnames(df)), drop = FALSE]
  cols_to_fix <- setdiff(colnames(df), c("pos", "start", "end"))
  df[cols_to_fix] <- lapply(df[cols_to_fix], function(x) {
    if (is.list(x) || is.factor(x)) as.character(x) else x
  })
  as.data.frame(df)
})

patient_parallel_display_data <- reactive({
  df <- patient_filtered_parallel_data()
  if (is.null(df) || nrow(df) == 0) return(data.frame(Message = "No parallel data found"))
  
  existing <- intersect(parallel_base_columns, colnames(df))
  others <- setdiff(colnames(df), existing)
  df <- df[, c(existing, others), drop = FALSE]
  
  cols_to_fix <- setdiff(colnames(df), c("pos", "start", "end"))
  df[cols_to_fix] <- lapply(df[cols_to_fix], function(x) {
    if (is.list(x) || is.factor(x)) as.character(x) else x
  })
  as.data.frame(df)
})