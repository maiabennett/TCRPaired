

getCDR3SeqCluster <- function(data, combine.chains = TRUE, filter.by = NULL) {
    if (combine.chains) {
        data <- data %>% 
            mutate(CDR3 = paste0(CDR3a, CDR3b)) %>% 
            getSeqCluster("CDR3", "AB", filter.by = filter.by)
    } else {
        data <- data %>% 
            getSeqCluster("CDR3a", "A", filter.by = filter.by) %>% 
            getSeqCluster("CDR3b", "B", filter.by = filter.by)
    }
    return(data)
}

getCDRSeqCluster <- function(data, combine.chains = TRUE, filter.by = NULL) {
    if (combine.chains) {
        data <- data %>%
            mutate(CDR.seq = paste0(CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b)) %>%
            getSeqCluster("CDR.seq", "AB", filter.by = filter.by)
    } else {
        data <- data %>% 
            mutate(CDRa.seq = paste0(CDR1a, CDR2a, CDR2.5a, CDR3a),
                   CDRb.seq = paste0(CDR1b, CDR2b, CDR2.5b, CDR3b)) %>%
            getSeqCluster("CDRa.seq", "A", filter.by = filter.by) %>%
            getSeqCluster("CDRb.seq", "B", filter.by = filter.by)
    }
    return(data)
}

getFullSeqCluster <- function(data, combine.chains = TRUE, filter.by = NULL) {
    data <- data %>%
        getSeqCluster("full.seq", "AB", filter.by = filter.by)
    return(data)
}

getSeqCluster <- function(data, column, chain, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    data <- data %>%
        filter(!str_detect(!!sym(column), "[^ACDEFGHIKLMNPQRSTVWY]")) %>%
        mutate(chain = chain) %>%
        djvdj::cluster_sequences(data_col = column, chain_col = "chain")
    data <- data %>%
        select(-chain)

    return(data)
}

