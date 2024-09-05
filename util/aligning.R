

use_condaenv("./util/.conda", required = TRUE)
source_python("./util/extractSequences.py")


###################################
# Sequence from alignment methods #
###################################

# Get full TCR sequences without writing to fasta
getFullSeq <- function(data) {
  data.in <- formatVJCDR3(data) %>%
    na.omit() %>%
    dplyr::filter_all(all_vars(. != ""))
  tcrs.full <- manage_sequences(data.in) %>% 
    dplyr::rename(alpha.seq = Alpha, 
           beta.seq = Beta) %>% 
    dplyr::mutate(full.seq = paste0(alpha.seq, beta.seq))
  tcrs.full.constant <- manage_sequences(data.in, constant_regions = TRUE) %>% 
    dplyr::rename(alpha.seq.constant = Alpha, 
           beta.seq.constant = Beta) %>%
    dplyr::mutate(full.seq.constant = paste0(alpha.seq.constant, beta.seq.constant))

  data <- bindSeqs(data, tcrs.full, tcrs.full.constant)
}

# Get CDR sequences
# If alleles = TRUE, append alleles ONLY for the purpose of CDR2.5 assignment
getCDRSeq <- function(data, alleles = FALSE) {
  data.in <- data %>%
    dplyr::select(clone.id, AV, BV) %>%
    na.omit() %>%
    dplyr::filter_all(all_vars(. != ""))
  if (alleles) {
    data.in <- data.in %>%
      dplyr::mutate(AV = paste0(AV, "*01"), BV = paste0(BV, "*01"))
  }
  tcrs.cdrs <- convert_data_paired(data.in, "AV", "BV")
  
  data <- bindSeqs(data, tcrs.cdrs)
}

# Bind new sequence entries to data while checking if they already exist (mostly for CDR sequences)
bindSeqs <- function(data, ...) {
  new.data <- list(...)
  for (new in new.data) {
    # Join the new data to the existing data
    data <- dplyr::left_join(data, new, by = "clone.id", suffix = c("", ".new"))
    # Get the shared column names
    common_cols <- intersect(names(data), paste0(names(new), ".new"))
    common_cols <- sub("\\.new$", "", common_cols)
    # Iterate over the shared columns
    for (col in common_cols) {
      # Replace NA or empty values in data with the corresponding values from new
      data <- data %>%
        dplyr::mutate(!!col := case_when(
          is.na(!!rlang::sym(col)) | !!rlang::sym(col) == "" ~ !!rlang::sym(paste0(col, ".new")),
          TRUE ~ !!rlang::sym(col)
        ))
    }
    # Remove the columns from new
    data <- data %>% dplyr::select(-ends_with(".new"))
  }
  return(data)
} 
