
getCDRCluster <- function(data, column, chain) {
    chain <- data %>%
        mutate(chain = chain) %>%
        cluster_sequences(data_col = column, chain_col = "chain")
    data <- chain %>%
        select(-chain)
    return(data)
}

