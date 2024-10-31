


######################
# Formatting methods #
######################
# Convert MHC and VDJ gene names to a standard format
convertGenes <- function(data) {

    gene.file <- "./util/ref/family_alignment.txt"
    gene.ref <- data.frame(ID = readChar(gene.file, file.info(gene.file)$size))
    gene.ref <- gene.ref %>% 
        dplyr::mutate(ID = strsplit(as.character(ID), ">")) %>%
        unnest(ID)
    gene.ref <- gene.ref[-1,]
    gene.ref <- gene.ref %>% 
        separate(ID, sep = "[|]", into=c("ID", "Name", "Sequence"), extra="merge", fill="right") %>%
        separate(Sequence, sep = "\\|(?=[^|]*$)", into=c("Rest", "Sequence"), extra="merge", fill="right", remove=FALSE) %>%
        select(Name, Sequence)
        
    genes.single <- gene.ref %>%
        dplyr::mutate(gene.name = str_extract(Name, "^[^\\*]*")) %>%
        # Group by Sequence to ensure unique sequences
        dplyr::group_by(gene.name) %>%
        dplyr::mutate(unique_seq_count = n_distinct(Sequence)) %>%
    # Filter to keep only those gene.names where there is exactly one unique sequence
        filter(unique_seq_count == 1) %>%
        dplyr::ungroup() %>%
        # Select distinct SimpleName-Sequence pairs
        dplyr::distinct(gene.name, Sequence, .keep_all = TRUE) %>%
        # Pull the SimpleName for the genes.single list
        dplyr::pull(gene.name)

    data <- data %>%
        # Remove murine MHC data
        filter(!str_detect(MHC, "H-2K")) %>%
        # fix gene naming (e.g., TCRAV --> TRAV)
        # make gene columns without alleles (e.g., TRAV14-1*01 --> TRAV14)
        # add allele only to single-allele genes (e.g., TRAV18 --> TRAV18*01)
        dplyr::mutate(AV = gsub("TCRA", "TRA", AV),
            AJ = gsub("TCRA", "TRA", AJ),
            BV = gsub("TCRB", "TRB", BV),
            BJ = gsub("TCRB", "TRB", BJ),
            AV.gene = str_extract(AV, "^[^\\*]*"),
            AJ.gene = str_extract(AJ, "^[^\\*]*"),
            BV.gene = str_extract(BV, "^[^\\*]*"),
            BJ.gene = str_extract(BJ, "^[^\\*]*"),
            AV.family = ifelse(str_detect(AV.gene, "/"),
                (str_replace_all(AV.gene, "(.*?)(/.*)", "\\1") %>%
                    str_extract("^[^\\-*]+")),
                str_extract(AV.gene, "^[^\\-*]+")),
            BV.family = ifelse(str_detect(BV.gene, "/"),
                (str_replace_all(BV.gene, "(.*?)(/.*)", "\\1") %>%
                    str_extract("^[^\\-*]+")),   
                str_extract(BV.gene, "^[^\\-*]+")),
            AJ.family = ifelse(str_detect(AJ.gene, "/"),
                (str_replace_all(AJ.gene, "(.*?)(/.*)", "\\1") %>%
                    str_extract("^[^\\-*]+")),
                str_extract(AJ.gene, "^[^\\-*]+")),
            BJ.family = ifelse(str_detect(BJ.gene, "/"),
                (str_replace_all(BJ.gene, "(.*?)(/.*)", "\\1") %>%
                    str_extract("^[^\\-*]+")),
                str_extract(BJ.gene, "^[^\\-*]+")),
            AV = ifelse(AV.gene %in% genes.single & !grepl("\\*", AV), paste0(AV, "*01"), AV),
            BV = ifelse(BV.gene %in% genes.single & !grepl("\\*", BV), paste0(BV, "*01"), BV),
            AJ = ifelse(AJ.gene %in% genes.single & !grepl("\\*", AJ), paste0(AJ, "*01"), AJ),
            BJ = ifelse(BJ.gene %in% genes.single & !grepl("\\*", BJ), paste0(BJ, "*01"), BJ),
            MHC = ifelse(!is.na(MHC) & MHC != "" & !startsWith(MHC, "H"), paste0("HLA-", MHC), MHC),
            MHC = str_remove(MHC, ",.*"),
            MHC = ifelse(str_detect(MHC, "\\*"), MHC, str_replace_all(MHC, "(\\D)(\\d)", "\\1*\\2")),
            MHC.allele = ifelse(str_detect(MHC, ":"), str_extract(MHC, "^[^:]+"), MHC),
            MHC.locus = str_extract(MHC, "HLA-[^\\d\\*w]+")
            )
    return(data)
}

# Convert epitope species names to a standard format
convertEpitopeSpecies <- function(data) {
    swap_values <- c("NY-ESO-1", "ORF1ab", "PIK3CA", "p53")
    data <- data %>% 
        dplyr::mutate(
        temp = ifelse(Epitope.species %in% swap_values, Epitope.species, NA),
        Epitope.species = ifelse(Epitope.species %in% swap_values, Epitope.gene, Epitope.species),
        Epitope.gene = ifelse(!is.na(temp), temp, Epitope.gene)
        ) %>%
        dplyr::select(-temp) 
        data <- data %>%
        dplyr::mutate(Epitope.species = case_when(str_detect(Epitope.species, "CMV|Human herpesvirus 5") ~ "Cytomegalovirus (CMV)",
        str_detect(Epitope.species, "EBV|Human herpesvirus 4") & !str_detect(Epitope.species, "Influenza") ~ "Epstein Barr virus (EBV)",
        str_detect(Epitope.species, "DENV|dengue") ~ "Dengue virus (DENV)",
        str_detect(Epitope.species, "H1N1") ~ "Influenza A (H1N1)",
        str_detect(Epitope.species, "H9N2") ~ "Influenza A (H9N2)",
        str_detect(Epitope.species, "HCV|Hepatitis C") ~ "Hepatitis C virus (HCV)",
        str_detect(Epitope.species, "HKU1") ~ "Human coronavirus HKU1 (HCoV-HKU1)",
        str_detect(Epitope.species, "OC43") ~ "Human coronavirus OC43 (HCoV-OC43)",
        str_detect(Epitope.species, "229E") ~ "Human coronavirus 229E (HCoV-229E)",
        str_detect(Epitope.species, "HIV-1") ~ "Human immunodeficiency virus 1 (HIV-1)",
        str_detect(Epitope.species, "HPV|papillomavirus") ~ "Human papillomavirus (HPV)",
        str_detect(Epitope.species, "HSV-2|HSV2|Human herpesvirus 2") ~ "Herpes simplex virus 2 (HSV-2)",
        str_detect(Epitope.species, "herpes virus 1|herpesvirus1") ~ "Herpes simplex virus 1 (HSV-1)",
        str_detect(Epitope.species, "HTLV-1|T cell leukemia virus") ~ "Human T-lymphotropic virus type 1 (HTLV-1)",
        str_detect(Epitope.species, "hepatitis B") ~ "Hepatitis B virus (HBV)",
        str_detect(Epitope.species, "Hepatitis E") ~ "Hepatitis E virus (HEV)",
        str_detect(Epitope.species, "Homo sapiens|HomoSapiens|Homo Sapiens|Neoantigen|Melanoma") ~ "Homo sapiens (human)",
        str_detect(Epitope.species, "Influenza A|InfluenzaA") & !str_detect(Epitope.species, "H1N1|H9N2|EBV") ~ "Influenza A",
        (Epitope.species == "Influenza") & str_detect(Epitope.gene, "M1") ~ "Influenza A",
        str_detect(Epitope.species, "Influenza B") ~ "Influenza B",
        str_detect(Epitope.species, "tuberculosis|Mtb") ~ "Mycobacterium tuberculosis",
        str_detect(Epitope.species, "MCPyV|Merkel cell") ~ "Merkel cell polyomavirus (MCPyV)",
        str_detect(Epitope.species, "MERS") ~ "Middle East respiratory syndrome coronavirus (MERS-CoV)",
        str_detect(Epitope.species, "BJ01") ~ "Severe acute respiratory syndrome coronavirus BJ01 (SARS-CoV-1 BJ01)",
        str_detect(Epitope.species, "Tor2") ~ "Severe acute respiratory syndrome coronavirus Tor2 (SARS-CoV-1 Tor2)",
        str_detect(Epitope.species, "Urbani") ~ "Severe acute respiratory syndrome coronavirus Urbani (SARS-CoV-1 Urbani)",
        str_detect(Epitope.species, "SARS-CoV-2|SARS-CoV2|coronavirus 2") ~ "Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)",
        str_detect(Epitope.species, "SARS-CoV1| coronavirus 1") & !str_detect(Epitope.species, "BJ01|Tor2|Urbani") ~ "Severe acute respiratory syndrome coronavirus 1 (SARS-CoV-1)",
        str_detect(Epitope.species, "YFV|Yellow fever") ~ "Yellow fever virus (YFV)",
        str_detect(Epitope.species, "StreptomycesKanamyceticus") ~ "Streptomyces kanamyceticus",
        str_detect(Epitope.species, "Vaccinia virus (vaccinia virus VV)") ~ "Vaccinia virus (VV)",
        str_detect(Epitope.species, "EBV") & str_detect(Epitope.species, "Influenza") ~ "Epstein Barr virus (EBV), Influenza A",
        TRUE ~ Epitope.species
        ))%>%
        dplyr::group_by(Epitope) %>%
        dplyr::mutate(Epitope.species = ifelse(Epitope.species == "" | is.na(Epitope.species), min(Epitope.species[nchar(Epitope.species) > 0], na.rm = TRUE), Epitope.species)) %>%
        dplyr::ungroup() %>%
        # Filter out epitope species erroneously in Epitope column
        filter(str_detect(Epitope, "^[A-Z]+$"))
}

# Convert epitope gene names to a standard format
convertEpitopeGenes <- function(data) {
    data <- data %>%
        dplyr::mutate(Epitope.gene = str_remove_all(Epitope.gene, "\\[.*?\\]")) %>%
        dplyr::mutate(Epitope.gene = case_when(str_detect(Epitope.gene, "BMLF|BMLF1") & !str_detect(Epitope.gene, "M1") ~ "BMLF1",
        str_detect(Epitope.gene, "BRLF1|BRLF-1") ~ "BRLF1",
        str_detect(Epitope.gene, "BZLF1|BZLF-1") ~ "BZLF1",
        str_detect(Epitope.gene, "M1|Matrix protein 1|matrix protein 1") & (Epitope.species == "Influenza A") ~ "M1",
        (Epitope.gene == "M") & str_detect(Epitope.species,"Influenza A") ~ "M1",
        str_detect(Epitope.gene, "BMLF") & str_detect(Epitope.gene, "M1") ~ "BMLF1, M1",
        str_detect(Epitope.gene, "EBNA1|EBNA-1|Epstein-Barr nuclear antigen 1")  ~ "EBNA1",
        str_detect(Epitope.gene, "EBNA3|EBNA-3") & !str_detect(Epitope.gene, "3A|3B|3C")  ~ "EBNA3",
        str_detect(Epitope.gene, "EBNA3A|EBNA-3A")  ~ "EBNA3A",
        str_detect(Epitope.gene, "EBNA3B|EBNA-3B")  ~ "EBNA3B",
        str_detect(Epitope.gene, "EBNA3C|EBNA-3C")  ~ "EBNA3C",
        str_detect(Epitope.gene, "EBNA4|EBNA-4")  ~ "EBNA4",
        str_detect(Epitope.gene, "EBNA6|EBNA-6")  ~ "EBNA6",
        str_detect(Epitope.gene, "GAG|Gag|gag") & !str_detect(Epitope.gene, "Pol |pol ")  ~ "Gag",
        str_detect(Epitope.gene, "NY-ESO-1")  ~ "NY-ESO-1",
        str_detect(Epitope.gene, "Nef")  ~ "Nef",
        str_detect(Epitope.gene, "Spike|spike|surface glycoprotein")  ~ "Spike",
        str_detect(Epitope.gene, "pol |Pol ") & !str_detect(Epitope.gene, "Gag|gag") ~ "Pol",
        str_detect(Epitope.gene, "LMP2A|LMP-2A|Latent membrane protein 2")  ~ "LMP2A",
        str_detect(Epitope.gene, "E6|Protein E6")  ~ "E6",
        str_detect(Epitope.gene, "E7|Protein E7")  ~ "E7",
        str_detect(Epitope.gene, "ORF1ab|orf1ab|Replicase polyprotein 1")  ~ "ORF1ab",
        str_detect(Epitope.gene, "ORF3|ORF3a|Replicase polyprotein 3")  ~ "ORF3a",
        str_detect(Epitope.gene, "ORF6|Replicase polyprotein 6")  ~ "ORF6",
        str_detect(Epitope.gene, "ORF7a|Replicase polyprotein 7a")  ~ "ORF7a",
        str_detect(Epitope.gene, "ORF7b|Replicase polyprotein 7b")  ~ "ORF7b",
        str_detect(Epitope.gene, "ORF8|Replicase polyprotein 8")  ~ "ORF8",
        str_detect(Epitope.gene, "ORF9b|Replicase polyprotein 9")  ~ "ORF9b",
        str_detect(Epitope.gene, "ORF10|Replicase polyprotein 10")  ~ "ORF10",
        str_detect(Epitope.gene, "Membrane protein|Matrix|M protein")  ~ "Membrane protein",
        TRUE ~ Epitope.gene)) %>%
        dplyr::group_by(Epitope) %>%
        dplyr::mutate(Epitope.gene = ifelse(Epitope.gene == "" | is.na(Epitope.gene), min(Epitope.gene[nchar(Epitope.gene) > 0]), Epitope.gene)) %>%
        dplyr::mutate(Epitope.gene = min(Epitope.gene[nchar(Epitope.gene) > 0])) %>%
        dplyr::ungroup()
}

# Identify epitopes above a specified number of receptors, and optionally join with epitope information from the reference dataset
fetchEpitopes <- function(data, column, threshold = 10, reference = TRUE) {
    epitopes <- data %>% 
        separate_rows(!!sym(column), sep = ",") %>%
        group_by(!!sym(column)) %>%
        count() %>%
        filter(n >= threshold)
    if (reference) {
        epitopes <- epitopes %>%
        left_join(data %>% 
            select(Epitope, Epitope.gene, Epitope.species) %>% 
            distinct(),
        by = column)
    }
    return(epitopes)
}

# Check if sequence is proper AA sequence (for CD-HIT alignment)
checkAASequence <- function(data) {
    aa <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", 
                        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
    data <- data %>% 
        filter(!str_detect(CDR3a, paste0("[^", paste(aa, collapse = ""), "]")), 
            !str_detect(CDR3b, paste0("[^", paste(aa, collapse = ""), "]")))
}

# Format full, paired dataset (i.e., start new benchmark compatible with all methods)
formatFullPaired <- function(data) {
    data <- data %>% 
        convertGenes() %>%
        convertEpitopeSpecies() %>%
        convertEpitopeGenes() %>%
        getFullSeq() %>%
        getCDRSeq() %>%
        filterFullSeq() %>%
        filterCDRSeq()
}

# Convert to V/J gene + CDR3 format (Train/SeqConductor)
formatVJCDR3 <- function(data) {
    #data <- filterVJCDR3(data)
    data <- data %>%
        select(clone.id, AV, CDR3a, AJ, BV, CDR3b, BJ)
}

# Convert to V/J gene + CDR3 + Epitope format (ERGOII)
    formatERGOII <- function(data, epitope) {
    data <- filterVJCDR3(data)
    data <- data %>%
        select(clone.id, AV, CDR3a, AJ, BV, CDR3b, BJ) %>% 
        rename(TRAV = AV, TRA = CDR3a, TRAJ = AJ, TRBV = BV, TRB = CDR3b, TRBJ = BJ) %>%
        mutate(Peptide = epitope,
            `T-Cell-Type` = "",
            MHC = "") %>%
        select(TRA, TRB, TRAV, TRAJ, TRBV, TRBJ, `T-Cell-Type`, Peptide, MHC, clone.id)
}

# Convert to SwarmTCR format
formatSwarmTCR <- function(data, flag = NULL, keep_epitopes = FALSE) {
    data <- filterCDRSeq(data)
    if (!is.null(flag)) {
        data <- data %>%
        mutate(FLAG = flag) %>%
        rename(TCR_ID = clone.id) %>% 
        filter(across(contains("CDR"), ~ !grepl("[^[:alnum:]]", .))) # remove rows with symbols in CDR sequences
        if (keep_epitopes) {
        data <- data %>%
        select(TCR_ID, FLAG, CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b, Epitope)
        } else {
        data <- data %>%
            select(TCR_ID, FLAG, CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b)
        }
    } else {
        data <- data %>%
        rename(TCR_ID = clone.id) %>% 
        filter(across(contains("CDR"), ~ !grepl("[^[:alnum:]]", .))) # remove rows with symbols in CDR sequences
        if (keep_epitopes) {
        data <- data %>%
            select(TCR_ID, CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b, Epitope)
        } else {
        data <- data %>%
            select(TCR_ID, CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b)
        }
    }
}

# Convert to TCRdist2 format
formatTCRdist2 <- function(data) {
    data <- data %>%
        mutate(subject = source.ID,
                id = clone.id,
                epitope = Epitope,
                a_nuqseq = alpha.seq,
                b_nuqseq = beta.seq,
                a_len = nchar(alpha.seq),
                b_len = nchar(beta.seq),
                a_quals = "",
                b_quals = "") 

    for (i in seq_len(nrow(data))){
        a_len <- data$a_len[i]
        b_len <- data$b_len[i]
        a_reps <- rep(50, a_len)
        b_reps <- rep(50, b_len)
        a_quals <- paste0(a_reps, collapse = ".")
        b_quals <- paste0(b_reps, collapse = ".")
        data$a_quals[i] <- a_quals
        data$b_quals[i] <- b_quals
    }

    data <- data %>%
        select(id, epitope, subject, a_nuqseq, b_nuqseq, a_quals, b_quals)
}