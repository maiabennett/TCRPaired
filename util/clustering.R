
getCDRCluster <- function(data, column, chain, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    data <- data %>%
        filter(!str_detect(!!sym(column), "[^ACDEFGHIKLMNPQRSTVWY]")) %>%
        mutate(chain = chain) %>%
        djvdj::cluster_sequences(data_col = column, chain_col = "chain")
    data <- chain %>%
        select(-chain)

    return(data)
}

