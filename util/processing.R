

######################
# Processing methods #
######################

# Basic import
import <- function(file) {
  if (grepl(".csv", file)) {
    data <- read.csv(file, encoding = "UTF-8")
  } else if (grepl(".txt", file)) {
    data <- read.delim(file, sep = " ")
    if (ncol(data) == 1) {
      data <- read.delim(file, sep = "\t")
    }
  } else if (grepl(".tsv", file)) {
    data <- read.delim(file, sep = "\t")
  } else if (grepl(".xlsx", file)) {
    data <- read_excel(file)
  } else {
    stop("File type not recognized")
  }
  return(data)
}

# Preprocess IEDB data
preprocessIEDB <- function(data) {
  data <- data %>%
    # select columns 
    select(Chain.1.4, Chain.1.3, 
           Chain.1.8, Chain.1.7,
           Chain.1.19, Chain.1.18,
           Chain.1.24, Chain.1.23,
           Chain.1.12, Chain.1.11,
           Chain.2.4, Chain.2.3, 
           Chain.2.8, Chain.2.7,
           Chain.2.19, Chain.2.18,
           Chain.2.24, Chain.2.23,
           Chain.2.12, Chain.2.11,
           Epitope.1, Epitope.2, Epitope.3,
           Assay.2,
           Assay.1, Reference) %>%
    # rename columns
    rename(AVcalc = Chain.1.4, AVcur = Chain.1.3, 
           AJcalc = Chain.1.8, AJcur = Chain.1.7,
           CDR1acalc = Chain.1.19, CDR1acur = Chain.1.18,
           CDR2acalc = Chain.1.24, CDR2acur = Chain.1.23,
           CDR3acalc = Chain.1.12, CDR3acur = Chain.1.11,
           BVcalc = Chain.2.4, BVcur = Chain.2.3, 
           BJcalc = Chain.2.8, BJcur = Chain.2.7,
           CDR1b = Chain.2.19, CDR1bcur = Chain.2.18,
           CDR2b = Chain.2.24, CDR2bcur = Chain.2.23,
           CDR3bcalc = Chain.2.12, CDR3bcur = Chain.2.11,
           Epitope = Epitope.1,
           Epitope.gene = Epitope.2, Epitope.species = Epitope.3,
           MHC = Assay.2,
           Source.ID = Assay.1) %>%
    # remove first row
    slice(-1) %>%
    # assign CDR3 sequences from calculated and curated columns
    mutate(CDR3a = case_when(CDR3acalc != "" ~ CDR3acalc,
                             CDR3acur != "" ~ CDR3acur,
                             TRUE ~ ""),
           CDR3b = case_when(CDR3bcalc != "" ~ CDR3bcalc,
                             CDR3bcur != "" ~ CDR3bcur,
                             TRUE ~ ""),
           AV = case_when(AVcur != "" ~ AVcur,
                          AVcalc != "" ~ AVcalc,
                          TRUE ~ ""),
           AJ = case_when(AJcur != "" ~ AJcur,
                          AJcalc != "" ~ AJcalc,
                          TRUE ~ ""),
           BV = case_when(BVcur != "" ~ BVcur,
                          BVcalc != "" ~ BVcalc,
                          TRUE ~ ""),
           BJ = case_when(BJcur != "" ~ BJcur,
                          BJcalc != "" ~ BJcalc,
                          TRUE ~ ""),
           CDR1a = case_when(CDR1acur != "" ~ CDR1acur,
                            CDR1acalc != "" ~ CDR1acalc,
                            TRUE ~ ""),
           CDR2a = case_when(CDR2acur != "" ~ CDR2acur,
                            CDR2acalc != "" ~ CDR2acalc,
                            TRUE ~ ""),
           CDR1b = case_when(CDR1bcur != "" ~ CDR1bcur,
                            CDR1b != "" ~ CDR1b,
                            TRUE ~ ""),
           CDR2b = case_when(CDR2bcur != "" ~ CDR2bcur,
                            CDR2b != "" ~ CDR2b,
                            TRUE ~ ""),
           # remove " + ..." from epitope
           Epitope = gsub(" \\+ .*", "", Epitope)) %>%
    # remove rows missing essential values
    filter_at(vars("AV", "CDR3a", "BV", "CDR3b", "Epitope"),
              all_vars(. != ""))
  data <- data %>% 
    # make clone.id
    group_by(Epitope) %>%
    mutate(clone.id = paste0("iedb_", Epitope, "_", row_number())) %>% 
    ungroup() %>%
    # add C...F wrapper for CDR3 sequence
    mutate(CDR3a = ifelse(grepl("^C.*F$", CDR3a), CDR3a, 
                          paste0("C", CDR3a, "F")),
           CDR3b = ifelse(grepl("^C.*F$", CDR3b), CDR3b, 
                          paste0("C", CDR3b, "F")),
           Score = "") %>% 
    select(clone.id, 
           AV, CDR1a, CDR2a, CDR3a, AJ,
           BV, CDR1b, CDR2b, CDR3b, BJ, 
           Epitope, Epitope.gene, Epitope.species,
           MHC, Source.ID, Reference, Score) %>% 
    distinct(AV, CDR3a, BV, CDR3b, .keep_all = TRUE)
}

# Preprocess VDJdb data
preprocessVDJ <- function(data) {
  data <- data %>%
      filter(complex.id != 0)
    a <- data %>%
      filter(Gene == 'TRA')
    b <- data %>%
      filter(Gene == 'TRB')
    data <- merge(a, b, by = 'complex.id', suffixes = c('.a', '.b')) %>%
      select('V.a', 'CDR3.a', 'J.a', 
             'V.b', 'CDR3.b', 'J.b', 
             'Epitope.a', 'Epitope.gene.a', 'Epitope.species.a', 
             'MHC.A.a', 'complex.id', 'Reference.a', 'Score.a') %>%
      rename(AV = V.a, CDR3a = CDR3.a, AJ = J.a,
             BV = V.b, CDR3b = CDR3.b, BJ = J.b,
             Epitope = Epitope.a, Epitope.gene = Epitope.gene.a, Epitope.species = Epitope.species.a,
             MHC = MHC.A.a, Source.ID = complex.id,
             Reference = Reference.a, Score = Score.a) %>%
      # remove rows missing essential values
      filter_at(vars("AV", "CDR3a", "BV", "CDR3b", "Epitope"),
                all_vars(. != ""))
    data <- data %>% 
      # make clone.id
      group_by(Epitope) %>%
      mutate(clone.id = paste0("vdj_", Epitope, "_", row_number())) %>% 
      ungroup() %>%
      # add empty columns for CDR1 + CDR2 sequences which are not provided by VDJdb
      mutate(CDR1a = "", CDR2a = "",
             CDR1b = "", CDR2b = "") %>% 
      select(clone.id,
             AV, CDR1a, CDR2a, CDR3a, AJ,
             BV, CDR1b, CDR2b, CDR3b, BJ, 
             Epitope, Epitope.gene, Epitope.species,
             MHC, Source.ID, Reference, Score) %>% 
      distinct(AV, CDR3a, BV, CDR3b, .keep_all = TRUE)
}

# Preprocess McPAS-TCR data
preprocessMcPAS <- function(data) {
  data <- data %>% 
    select('TRAV', 'CDR3.alpha.aa', 'TRAJ',
           'TRBV', 'CDR3.beta.aa', 'TRBJ',
           'Epitope.peptide', 'Antigen.protein', 'Pathology',
           'MHC', 'PubMed.ID', 'Species') %>%
    rename(AV = TRAV, CDR3a = `CDR3.alpha.aa`, AJ = TRAJ,
           BV = TRBV, CDR3b = `CDR3.beta.aa`, BJ = TRBJ,
           Epitope = `Epitope.peptide`, Epitope.gene = Antigen.protein,
           Epitope.species = Pathology,
           Reference = PubMed.ID) %>%
    filter(Species == "Human") %>%
    # remove rows missing essential values
    filter_at(vars("AV", "CDR3a", "BV", "CDR3b", "Epitope"),
              all_vars(. != ""))
  data <- data %>%
    # make clone.id
    group_by(Epitope) %>%
    mutate(clone.id = paste0("mcpas_", Epitope, "_", row_number())) %>% 
    ungroup() %>%
    # add empty columns for CDR1 + CDR2 sequences which are not provided by McPAS-TCR
    # also, convert HLA naming to IMGT-compatible format 
    mutate(across(everything(), ~ iconv(., from = "UTF-8", to = "ASCII//TRANSLIT"))) %>%
    mutate(CDR1a = "", CDR2a = "",
           CDR1b = "", CDR2b = "",
           Score = "", Source.ID = "",
           #AV = iconv(AV, from = "ISO-8859-1", to = "UTF-8"),
           #MHC = iconv(MHC, from = "ISO-8859-1", to = "UTF-8"),
           AV = sub(' .*', '', gsub(':', '*', AV)),
           AJ = sub(' .*', '', gsub(':', '*', AJ)),
           BV = sub(' .*', '', gsub(':', '*', BV)),
           BJ = sub(' .*', '', gsub(':', '*', BJ)),
           #CDR3a = sub('l', 'L', CDR3a),
           #CDR3a = sub("�", "", CDR3a),
           #CDR3b = sub('l', 'L', CDR3b),
           #CDR3a = sub("�", "", CDR3b)
           ) %>%
    select(clone.id,
           AV, CDR1a, CDR2a, CDR3a, AJ,
           BV, CDR1b, CDR2b, CDR3b, BJ, 
           Epitope, Epitope.gene, Epitope.species,
           MHC, Source.ID, Reference, Score) %>% 
    distinct(AV, CDR3a, BV, CDR3b, .keep_all = TRUE)
}

# Where epitope indicates if there is an epitope and name is the name to use in clone.id
# Currently, assumes no epitope
preprocess10X <- function(data, epitope, name) {
  name <- ifelse(is.null(name), "10x", name)
  counts <- data %>%
    count(barcode, sort=TRUE, name='count') %>% 
    filter(count == 2)
  barcodes <- as.list(counts$barcode)
  data <- data %>%
    filter(barcode %in% barcodes)

  a <- data[data$chain == 'TRA', ]
  b <- data[data$chain == 'TRB', ]
  data <- merge(a, b, by = 'barcode', suffixes = c('.a', '.b')) 

  data <- data %>% 
    mutate(clone.id = paste0(name, "_", row_number())) %>%
    select(clone.id, barcode, 
      v_gene.a, cdr1.a, cdr2.a, cdr3.a, j_gene.a,
      v_gene.b, cdr1.b, cdr2.b, cdr3.b, j_gene.b,
      fwr1.a, fwr2.a, fwr3.a, fwr4.a,
      fwr1.b, fwr2.b, fwr3.b, fwr4.b,
      fwr1_nt.a, cdr1_nt.a, fwr2_nt.a, cdr2_nt.a, fwr3_nt.a, cdr3_nt.a, fwr4_nt.a,
      fwr1_nt.b,  cdr1_nt.b, fwr2_nt.b, cdr2_nt.b, fwr3_nt.b, cdr3_nt.b, fwr4_nt.b,
      raw_clonotype_id.a, reads.a, umis.a, reads.b, umis.b,
      length.a, length.b) %>%
    rename(AV = v_gene.a, CDR1a = cdr1.a, CDR2a = cdr2.a, CDR3a = cdr3.a, AJ = j_gene.a,
           BV = v_gene.b, CDR1b = cdr1.b, CDR2b = cdr2.b, CDR3b = cdr3.b, BJ = j_gene.b,
           FWR1a = fwr1.a, FWR2a = fwr2.a, FWR3a = fwr3.a, FWR4a = fwr4.a,
           FWR1b = fwr1.b, FWR2b = fwr2.b, FWR3b = fwr3.b, FWR4b = fwr4.b,
           FWR1a_nt = fwr1_nt.a, CDR1a_nt = cdr1_nt.a, FWR2a_nt = fwr2_nt.a, CDR2a_nt = cdr2_nt.a, FWR3a_nt = fwr3_nt.a, CDR3a_nt = cdr3_nt.a, FWR4a_nt = fwr4_nt.a,
           FWR1b_nt = fwr1_nt.b, CDR1b_nt = cdr1_nt.b, FWR2b_nt = fwr2_nt.b, CDR2b_nt = cdr2_nt.b, FWR3b_nt = fwr3_nt.b, CDR3b_nt = cdr3_nt.b, FWR4b_nt = fwr4_nt.b,
           clonotype = raw_clonotype_id.a) %>%
    mutate(alpha.seq = paste0(FWR1a, CDR1a, FWR2a, CDR2a, FWR3a, CDR3a, FWR4a),
           beta.seq = paste0(FWR1b, CDR1b, FWR2b, CDR2b, FWR3b, CDR3b, FWR4b),
           alpha.nucseq = paste0(FWR1a_nt, CDR1a_nt, FWR2a_nt, CDR2a_nt, FWR3a_nt, CDR3a_nt, FWR4a_nt),
           beta.nucseq = paste0(FWR1b_nt, CDR1b_nt, FWR2b_nt, CDR2b_nt, FWR3b_nt, CDR3b_nt, FWR4b_nt))
  

}

# Preprocess data based on type
preprocess <- function(data, type, epitope = FALSE, name = NULL) {
  if (type == "IEDB") {
    preprocessIEDB(data)
  } 
  else if (type == "VDJ") {
    preprocessVDJ(data)
  }
  else if (type == "McPAS") {
    preprocessMcPAS(data)
  }
  else if (type == "10X") {
    preprocess10X(data, epitope, name)
  }
  else {
    stop("Type not recognized")
  }
}

