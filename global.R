# global.R
library(shiny)
library(DT)
library(readr)
library(dplyr)
library(igvShiny)
library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)
library(GenomicFeatures) 
library(data.table)

if (requireNamespace("txdbmaker", quietly = TRUE)) {
  library(txdbmaker)
}

# ----------------- 1. FILE SYSTEM HELPERS -----------------
addResourcePath("splicing_plots", file.path(getwd(), "R_FRASER"))

find_file_robust <- function(relative_path) {
  if (file.exists(relative_path)) return(relative_path)
  p2 <- file.path("app", relative_path)
  if (file.exists(p2)) return(p2)
  p3 <- file.path("..", relative_path)
  if (file.exists(p3)) return(p3)
  return(NULL)
}

load_data_robust <- function(relative_path, load_func, ...) {
  found_path <- find_file_robust(relative_path)
  if (!is.null(found_path)) {
    message("Loading: ", found_path)
    tryCatch({
      return(load_func(found_path, ...))
    }, error = function(e) {
      message("ERROR loading ", relative_path, ": ", e$message)
      return(NULL)
    })
  } else {
    warning("FILE NOT FOUND: ", relative_path)
    return(NULL)
  }
}

# ----------------- 2. BAM FILE INFRASTRUCTURE -----------------

# 1. Create Temp Directory for Slices
temp_bam_dir <- file.path("www", "tmp_bam_slices")
if (!dir.exists(temp_bam_dir)) dir.create(temp_bam_dir, recursive = TRUE)

# 2. SMART RECURSIVE BAM FINDER
find_bam_file <- function(subfolder_hint, filename) {
  # A. Try Direct Hints (Standard & WWW)
  possible_paths <- c(
    file.path(subfolder_hint, filename),
    file.path("Input", subfolder_hint, filename),
    file.path("www", subfolder_hint, filename),      
    file.path("app", "www", subfolder_hint, filename) 
  )
  
  for (path in possible_paths) {
    if (file.exists(path)) {
      # Check index
      if (file.exists(paste0(path, ".bai")) || file.exists(sub(".bam$", ".bai", path))) {
        message("BAM Found (Direct): ", path)
        return(path)
      }
    }
  }
  
  # B. Recursive Search (Deep Scan)
  message("Searching recursively for BAM: ", filename, " ...")
  
  # Define search roots
  roots <- c(".", "Input", "www", "../Input", "app/Input", "app/www")
  roots <- roots[dir.exists(roots)]
  
  for (root in roots) {
    hits <- list.files(root, pattern = paste0("^", filename, "$"), recursive = TRUE, full.names = TRUE)
    
    for (hit in hits) {
      bai_1 <- paste0(hit, ".bai")
      bai_2 <- sub(".bam$", ".bai", hit)
      
      if (file.exists(bai_1) || file.exists(bai_2)) {
        message("BAM Found (Recursive): ", hit)
        return(hit)
      }
    }
  }
  
  warning("BAM NOT FOUND (or missing .bai index): ", filename)
  return(NULL)
}

# ----------------- 3. DATA LOADING -----------------
rna_path <- find_file_robust("Input/sample_pairs.txt")
if(!is.null(rna_path)) {
  rna_to_dna <- read.table(rna_path, col.names = c("rna_id", "dna_id"))
} else {
  rna_to_dna <- NULL; warning("Sample pairs file not found!")
}

var_path <- find_file_robust("FIN_TAB/annotated_variants.rds")
variants <- if(!is.null(var_path)) read_rds(var_path) else NULL

par_path <- find_file_robust("FIN_TAB/parallel_dna_rna.rds")
if(!is.null(par_path)) {
  parallel_data <- read_rds(par_path) %>% ungroup() %>%
    dplyr::select(any_of(c("Gene", "SYMBOL", "Feature", "HGVSc", "HGVSp",
                           "location", "Consequence", "chr", "pos", "strand",
                           "ref", "alt", "sample_id", "PARALLEL_STATUS", "DNA_ONLY_count",
                           "RNA_ONLY_count", "BOTH_count",
                           "rna_qual", "rna_GT", "rna_REF_AD", "rna_ALT_AD",
                           "rna_DP", "rna_GQ", "rna_PL", "rna_PGT", "rna_PID", "rna_PS",
                           "rna_VAF", "RNA_ONLY_count", "dna_qual", "dna_GT", "dna_REF_AD", "dna_ALT_AD",
                           "dna_DP", "dna_GQ", "dna_PL", "dna_PGT", "dna_PID", "dna_PS", "dna_VAF"
    ))) %>% as.data.frame()
} else { parallel_data <- NULL }

ase_path <- find_file_robust("R_ASE/ase_stat.csv")
if(!is.null(ase_path)) {
  ase_data <- readr::read_csv(ase_path, col_types = cols(.default = col_character())) %>%
    mutate(SYMBOL = as.character(SYMBOL)) %>% as.data.frame()
} else { ase_data <- NULL }

# 🚀 HIGH SPEED LOADING FOR OUTLIER TABLES
batch_counts <- load_data_robust("R_FEATURE_COUNTS/batch_counts.txt", data.table::fread, data.table = FALSE)
outrider_data <- load_data_robust("R_OUTRIDER/outrider.csv", data.table::fread, data.table = FALSE)

# --- FRASER LOADING ---
fraser_data <- load_data_robust("R_FRASER/fraser.csv", data.table::fread, data.table = FALSE)
if (is.null(fraser_data)) {
  fraser_data <- NULL
} else {
  colnames(fraser_data) <- tolower(colnames(fraser_data))
  cols <- colnames(fraser_data)
  chr_idx <- grep("^(chr|seqnames|chromosome|seqname)$", cols)[1]
  if (!is.na(chr_idx)) colnames(fraser_data)[chr_idx] <- "chr"
  start_idx <- grep("^(start|st)$", cols)[1]
  if (!is.na(start_idx)) colnames(fraser_data)[start_idx] <- "start"
  end_idx <- grep("^(end|en)$", cols)[1]
  if (!is.na(end_idx)) colnames(fraser_data)[end_idx] <- "end"
  
  as_num_safe <- function(x) { if(is.factor(x)) as.numeric(as.character(x)) else as.numeric(x) }
  if("start" %in% colnames(fraser_data)) fraser_data$start <- as_num_safe(fraser_data$start)
  if("end" %in% colnames(fraser_data))   fraser_data$end   <- as_num_safe(fraser_data$end)
  
  message("Fraser loaded. Columns: ", paste(colnames(fraser_data), collapse=", "))
}

# =====================================================================
# --- GLOBAL CLEANUP: Strip ".bam" from all patient/sample IDs ---
# =====================================================================

# 1. Parallel Data
if (!is.null(parallel_data) && "sample_id" %in% colnames(parallel_data)) {
  parallel_data$sample_id <- gsub("\\.bam$", "", parallel_data$sample_id, ignore.case = TRUE)
}

# 2. ASE Data
if (!is.null(ase_data)) {
  samp_col <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(ase_data), value = TRUE, ignore.case = TRUE)[1]
  if (!is.na(samp_col)) ase_data[[samp_col]] <- gsub("\\.bam$", "", ase_data[[samp_col]], ignore.case = TRUE)
}

# 3. Outrider Data
if (!is.null(outrider_data)) {
  samp_col <- grep("^(sample_id|sampleID|sample|rna_id|patient)$", colnames(outrider_data), value = TRUE, ignore.case = TRUE)[1]
  if (!is.na(samp_col)) outrider_data[[samp_col]] <- gsub("\\.bam$", "", outrider_data[[samp_col]], ignore.case = TRUE)
}

# 4. Fraser Data
if (!is.null(fraser_data)) {
  samp_col <- grep("^(sampleID|sample_id|sample|rna_id|patient)$", colnames(fraser_data), value = TRUE, ignore.case = TRUE)[1]
  if (!is.na(samp_col)) fraser_data[[samp_col]] <- gsub("\\.bam$", "", fraser_data[[samp_col]], ignore.case = TRUE)
}

# 5. Counts Data (Columns are sample IDs)
if (!is.null(batch_counts)) {
  colnames(batch_counts) <- gsub("\\.bam$", "", colnames(batch_counts), ignore.case = TRUE)
}

# 6. RNA-to-DNA mapping
if (!is.null(rna_to_dna)) {
  if ("rna_id" %in% colnames(rna_to_dna)) rna_to_dna$rna_id <- gsub("\\.bam$", "", rna_to_dna$rna_id, ignore.case = TRUE)
  if ("dna_id" %in% colnames(rna_to_dna)) rna_to_dna$dna_id <- gsub("\\.bam$", "", rna_to_dna$dna_id, ignore.case = TRUE)
}

# ----------------- 4. ANNOTATION CACHING -----------------

exons_gr <- NULL
introns_gr <- NULL

input_dir_path <- find_file_robust("Input")
if(is.null(input_dir_path)) input_dir_path <- "." 

cache_exons_path   <- file.path(input_dir_path, "exons_gr.rds")
cache_introns_path <- file.path(input_dir_path, "introns_gr.rds")

if (file.exists(cache_exons_path) && file.exists(cache_introns_path)) {
  message("Loading cached Annotation...")
  exons_gr <- readRDS(cache_exons_path)
  introns_gr <- readRDS(cache_introns_path)
  
} else {
  message("Cache not found. Loading GTF/TxDb (Slow ONCE)...")
  target_gtf_name <- "Homo_sapiens.GRCh38.114.gtf.gz"
  gtf_path <- find_file_robust(file.path("Input", target_gtf_name))
  
  if (is.null(gtf_path)) {
    candidates <- list.files(input_dir_path, pattern = "\\.gtf(\\.gz)?$", full.names = TRUE)
    if (length(candidates) > 0) gtf_path <- candidates[1]
  }
  
  txdb_path <- find_file_robust("Input/txdb.sqlite")
  save_txdb_path <- if (!is.null(gtf_path)) file.path(dirname(gtf_path), "txdb.sqlite") else "txdb.sqlite"
  
  txdb <- NULL
  if (!is.null(txdb_path)) {
    tryCatch({ txdb <- loadDb(txdb_path) }, error = function(e) { message("Error loading TxDb: ", e$message) })
  } 
  
  if (is.null(txdb) && !is.null(gtf_path)) {
    message("Parsing GTF: ", gtf_path)
    tryCatch({
      if (requireNamespace("txdbmaker", quietly = TRUE)) {
        txdb <- txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf")
      } else {
        txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_path, format = "gtf")
      }
      saveDb(txdb, save_txdb_path)
    }, error = function(e) { message("Error parsing GTF: ", e$message) })
  }
  
  if (!is.null(txdb)) {
    exons_gr <- exons(txdb, columns=c("gene_id"))
    introns_gr <- unlist(intronsByTranscript(txdb, use.names=TRUE))
    
    tryCatch({
      saveRDS(exons_gr, cache_exons_path)
      saveRDS(introns_gr, cache_introns_path)
      message("Cache saved!")
    }, error = function(e) { warning("Could not save cache: ", e$message) })
  } else {
    warning("CRITICAL: Could not create Annotation.")
  }
}

# ----------------- 5. CLEANUP TYPES -----------------
if(!is.null(variants)) {
  variants$chr <- as.character(variants$chr)
  variants$ref <- as.character(variants$ref)
  variants$alt <- as.character(variants$alt)
}
if(!is.null(parallel_data)) {
  parallel_data$chr <- as.character(parallel_data$chr)
  parallel_data$ref <- as.character(parallel_data$ref)
  parallel_data$alt <- as.character(parallel_data$alt)
}

genome_options <- parseAndValidateGenomeSpec(genomeName = "hg38", initialLocus = "BRCA1")

required_columns <- c("Gene", "SYMBOL", "Feature", "HGVSc", "HGVSp", "location", "Consequence", "chr", "pos", "strand", "ref", "alt", "AF", "CLIN_SIG", "Genomic_Context")
group_list <- list(
  "Variant Consequence" = c("Consequence", "IMPACT", "Genomic_Context"),
  "Sequence Context and Coding Potential" = c("cDNA_position", "CDS_position", "Protein_position",
                                              "Amino_acids", "Codons", "Existing_variation",
                                              "DISTANCE", "HGVSc", "HGVSp", "HGVS_OFFSET", "HGVS_IntronEndOffset", "HGVS_IntronStartOffset"), 
  "Gene Info" = c("chr", "pos", "SYMBOL", "location", "strand", "ref", "alt", "Gene",
                  "Feature", "Feature_type", "BIOTYPE", "CANONICAL", "MANE", "MANE_SELECT", "MANE_PLUS_CLINICAL", "RefSeq", "EXON", "INTRON"),
  "Frequencies" = c("AF", "gnomADg_AF", "gnomADg_AFR_AF", "gnomADg_AMI_AF", "gnomADg_AMR_AF",
                    "gnomADg_ASJ_AF", "gnomADg_EAS_AF", "gnomADg_FIN_AF", "gnomADg_MID_AF", "gnomADg_NFE_AF", "gnomADg_REMAINING_AF", "gnomADg_SAS_AF", "MAX_AF", "MAX_AF_POPS"),
  "Clinical Annotations" = c("ID", "GENE_PHENO", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED", "ClinVar", "ClinVar_ALLELEID", "ClinVar_CLNSIG", "ClinVar_CLNREVSTAT", "ClinVar_CLNDN", "ClinVar_CLNSIGCONF"),
  "RNA Metabolism" = c("SpliceAI_pred", "SpliceRegion", "SpliceVault_SpliceAI_delta", "SpliceVault_out_of_frame_events",
                       "SpliceVault_site_max_depth", "SpliceVault_site_pos", "SpliceVault_site_sample_count", "SpliceVault_site_type", "SpliceVault_top_events", "NMD", "NearestExonJB"),
  "Protein Impact" = c("SIFT", "PolyPhen", "DOMAINS", "am_class", "am_pathogenicity", "BayesDel", "CADD_PHRED", "CADD_RAW", "REVEL", "DownstreamProtein", "ProteinLengthChange"),
  "Functional Effect" = c("Enformer_SAD", "Enformer_SAR", "MaveDB_nt", "MaveDB_pro", "MaveDB_score", "MaveDB_urn", "pHaplo", "pTriplo"),
  "Regulatory Elements" = c("MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "TRANSCRIPTION_FACTORS",
                            "5UTR_annotation", "5UTR_consequence", "Existing_InFrame_oORFs", "Existing_OutOfFrame_oORFs", "Existing_uORFs"),
  "Technical Info" = c("FLAGS", "SYMBOL_SOURCE", "HGNC_ID", "REFSEQ_MATCH", "SOURCE", "REFSEQ_OFFSET", "GIVEN_REF", "USED_REF", "BAM_EDIT")
  
  
)

source("ui_helpers.R")
source("variants_module.R")
source("patient_module.R")

# =====================================================================
# --- 6. GLOBAL COLUMN DISPLAY & RENAMING CONFIGURATION ---
# =====================================================================

# 1. Rename Dictionary: Maps Database names to UI/Download names
# NOTE: The left side must EXACTLY match the database spelling!
column_names_mapping <- c(
  "gene" = "Gene",
  "symbol" = "SYMBOL"

)

# 2. Main Table Fixed Left-Side Columns
# The rest of the Main Table is dynamically ordered by the UI Checkboxes (which default to 'required_columns')
main_anchor_columns <- c("Gene", "SYMBOL", "Feature", "HGVSc", "HGVSp",
                         "location", "Consequence", "chr", "pos", "strand",
                         "ref", "alt")

# 4. Base order for PARALLEL Tables (Independent from Main Table)
parallel_base_columns <- c("Gene", "SYMBOL", "Feature", "HGVSc", "HGVSp",
                           "location", "Consequence", "chr", "pos", "strand",
                           "ref", "alt", "sample_id", "PARALLEL_STATUS", "DNA_ONLY_count",
                           "RNA_ONLY_count", "BOTH_count",
                           "rna_qual", "rna_GT", "rna_REF_AD", "rna_ALT_AD",
                           "rna_DP", "rna_GQ", "rna_PL", "rna_PGT", "rna_PID", "rna_PS",
                           "rna_VAF", "RNA_ONLY_count", "dna_qual", "dna_GT", "dna_REF_AD", "dna_ALT_AD",
                           "dna_DP", "dna_GQ", "dna_PL", "dna_PGT", "dna_PID", "dna_PS", "dna_VAF"
)

# 5. Cosmetic Helper Function
format_for_display <- function(df, sort_cols = NULL) {
  if (is.null(df) || nrow(df) == 0) return(df)
  
  # A. Rename FIRST (Applies the cosmetic UI names so we can sort by them!)
  cur_names <- colnames(df)
  new_names <- sapply(cur_names, function(n) {
    if (n %in% names(column_names_mapping)) column_names_mapping[[n]] else n
  })
  colnames(df) <- unname(new_names)
  
  # B. Reorder LAST (Matches against the beautiful UI names, Case-Insensitive!)
  if (!is.null(sort_cols)) {
    existing <- character(0)
    
    for (sc in sort_cols) {
      # Forgiving match: Ignores uppercase/lowercase differences
      match_idx <- which(tolower(colnames(df)) == tolower(sc))
      if (length(match_idx) > 0) {
        existing <- c(existing, colnames(df)[match_idx[1]])
      }
    }
    
    existing <- unique(existing)
    others <- setdiff(colnames(df), existing)
    df <- df[, c(existing, others), drop = FALSE]
  }
  
  return(df)
}

# 3. Base order for Outlier Tables
# (Any columns not listed here are just pushed to the far right automatically)
outlier_base_columns <- c("sample_id", "Gene", "SYMBOL")
ase_base_columns     <- c("sample_id", "Gene", "SYMBOL")
counts_base_columns  <- c("sample_id", "Gene", "SYMBOL")
fraser_base_columns  <- c("sample_id", "Gene", "SYMBOL")
