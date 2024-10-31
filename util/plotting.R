

## Plotting functions for V/J genes
plotVJAllele <- function(data, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    data <- data %>%
        select(clone.id, AV, AJ, BV, BJ, AV.gene, BV.gene, AJ.gene, BJ.gene, Epitope, MHC, Epitope.species) %>%
        pivot_longer(cols = c(AV, AJ, BV, BJ), names_to = "Gene Segment", values_to = "Allele") %>%
        mutate(keep = case_when(
            `Gene Segment` == "AV" ~ "AV.gene",
            `Gene Segment` == "AJ" ~ "AJ.gene",
            `Gene Segment` == "BV" ~ "BV.gene",
            `Gene Segment` == "BJ" ~ "BJ.gene"
        )) %>%
        pivot_longer(cols = c(AV.gene, AJ.gene, BV.gene, BJ.gene), names_to = "Gene Fam. Segment", values_to = "Gene") %>%
        filter(keep == `Gene Fam. Segment`) %>%
        group_by(Gene) %>%
        add_count(Allele, name = "Count") %>% 
        ungroup()

    plot <- ggplot(data,  
        aes(x = Allele, y = Count, fill = Gene)) +
        geom_bar(stat = "identity") +
        facet_wrap(~`Gene`, strip.position = "bottom", scales = "free", nrow = 1) +
        facet_wrap(~`Gene Segment`, scales = "free") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              strip.text = element_text(size = 8),
              legend.position = "none") +
        labs(title = "V/J Alleles",
                x = "Allele",
                y = "Count") +
        scale_fill_viridis_d(option = "viridis")

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}


plotVJGene <- function(data, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    data <- data %>%
            pivot_longer(cols = c(AV, AJ, BV, BJ), names_to = "Gene Segment", values_to = "Allele") %>%
            mutate(keep = case_when(
                `Gene Segment` == "AV" ~ "AV.gene",
                `Gene Segment` == "AJ" ~ "AJ.gene",
                `Gene Segment` == "BV" ~ "BV.gene",
                `Gene Segment` == "BJ" ~ "BJ.gene"
            )) %>%
            pivot_longer(cols = c(AV.gene, AJ.gene, BV.gene, BJ.gene), names_to = "Gene Fam. Segment", values_to = "Gene") %>%
            filter(keep == `Gene Fam. Segment`) %>%
            mutate(keep = case_when(
                `Gene Fam. Segment` == "AV.gene" ~ "AV.family",
                `Gene Fam. Segment` == "AJ.gene" ~ "AJ.family",
                `Gene Fam. Segment` == "BV.gene" ~ "BV.family",
                `Gene Fam. Segment` == "BJ.gene" ~ "BJ.family"
            )) %>%
            pivot_longer(cols = c(AV.family, AJ.family, BV.family, BJ.family), names_to = "Gene Chain", values_to = "Gene Family") %>%
            filter(keep == `Gene Chain`) %>%
            group_by(`Gene Family`) %>%
            add_count(Gene, name = "Count") %>% 
            ungroup()
    plot <- ggplot(data, 
        aes(x = Gene, fill = `Gene Family`)) +
        geom_bar(position = "dodge") +
        facet_wrap(~`Gene Family`, strip.position = "bottom", scales = "free", nrow = 1) +
        facet_wrap(~`Gene Segment`, scales = "free") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              strip.text = element_text(size = 8),
              legend.position = "none") +
        labs(title = "V/J Genes",
                x = "Gene",
                y = "Count") +
        scale_fill_viridis_d(option = "viridis")

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

plotVJFamily <- function(data, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    data <- data %>%
            pivot_longer(cols = c(AV, AJ, BV, BJ), names_to = "Gene Segment", values_to = "Allele") %>%
            mutate(keep = case_when(
                `Gene Segment` == "AV" ~ "AV.gene",
                `Gene Segment` == "AJ" ~ "AJ.gene",
                `Gene Segment` == "BV" ~ "BV.gene",
                `Gene Segment` == "BJ" ~ "BJ.gene"
            )) %>%
            pivot_longer(cols = c(AV.gene, AJ.gene, BV.gene, BJ.gene), names_to = "Gene Fam. Segment", values_to = "Gene") %>%
            filter(keep == `Gene Fam. Segment`) %>%
            mutate(keep = case_when(
                `Gene Fam. Segment` == "AV.gene" ~ "AV.family",
                `Gene Fam. Segment` == "AJ.gene" ~ "AJ.family",
                `Gene Fam. Segment` == "BV.gene" ~ "BV.family",
                `Gene Fam. Segment` == "BJ.gene" ~ "BJ.family"
            )) %>%
            pivot_longer(cols = c(AV.family, AJ.family, BV.family, BJ.family), names_to = "Gene Chain", values_to = "Gene Family") %>%
            filter(keep == `Gene Chain`) %>% 
            add_count(`Gene Family`, name = "Count") %>%
            ungroup()
    ggplot(data, 
        aes(x = `Gene Family`, fill = `Gene Family`)) +
        geom_bar(position = "dodge") +
        facet_wrap(~`Gene Family`, strip.position = "bottom", scales = "free", nrow = 1) +
        facet_wrap(~`Gene Segment`, scales = "free") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              strip.text = element_text(size = 8),
              legend.position = "none") +
        labs(title = "V/J Genes",
                x = "Gene",
                y = "Count") +
        scale_fill_viridis_d(option = "viridis")

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}


## Plotting functions for epitope
plotEpitopeCountDistribution <- function(data, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    plot <- plotCountDistribution(data, "Epitope") 

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

plotEpitopeCountDistributionByFactor <- function(data, factor, top = NULL, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }
    
    if (is.null(top)) {
        plot <- plotCountDistributionByFactor(data, "Epitope", factor)
    } else {
        plot <- plotTopCountDistributionByFactor(data, "Epitope", factor, top)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

plotEpitopeComposition <- function(data, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }
            
    plot <- plotComposition(data, "Epitope")

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

plotEpitopeCompositionByFactor <- function(data, factor, top = NULL, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }
    
    if (is.null(top)) {
        plot <- plotCompositionByFactor(data, "Epitope", factor)
    } else {
        plot <- plotTopCountDistributionByFactor(data, "Epitope", factor, top)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

plotEpitopeSpeciesComposition <- function(data, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    plot <- plotComposition(data, "Epitope.species")

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

plotEpitopeSpeciesCompositionByFactor <- function(data, factor, top = NULL, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }
    
    if (is.null(top)) {
        plot <- plotCompositionByFactor(data, "Epitope.species", factor)
    } else {
        plot <- plotTopCountDistributionByFactor(data, "Epitope.species", factor, top)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

## Plotting functions for MHC
plotMHCAlleleCountDistribution <- function(data, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    plot <- plotCountDistribution(data, "MHC") 

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

plotMHCAlleleCountDistributionByFactor <- function(data, factor, top = NULL, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }
    
    if (is.null(top)) {
        plot <- plotCountDistributionByFactor(data, "MHC", factor)
    } else {
        plot <- plotTopCountDistributionByFactor(data, "MHC", factor, top)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

plotMHCComposition <- function(data, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }
    
    plot <- plotComposition(data, "MHC")

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}


plotMHCCompositionByFactor <- function(data, factor, top = NULL, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }
    
    if (is.null(top)) {
        plot <- plotCompositionByFactor(data, "MHC", factor)
    } else {
        plot <- plotTopCountDistributionByFactor(data, "MHC", factor, top)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

plotMHCLocusComposition <- function(data, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    plot <- plotComposition(data, "MHC.locus")

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

plotMHCLocusCompositionByFactor <- function(data, factor, top = NULL, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }
    
    if (is.null(top)) {
        plot <- plotCompositionByFactor(data, "MHC.locus", factor)
    } else {
        plot <- plotTopCountDistributionByFactor(data, "MHC.locus", factor, top)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

plotMHCAlleleComposition <- function(data, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    plot <- plotComposition(data, "MHC.allele")

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

plotMHCAlleleCompositionByFactor <- function(data, factor, top = NULL, plotly = TRUE, filter.by = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }
    
    if (is.null(top)) {
        plot <- plotCompositionByFactor(data, "MHC.allele", factor)
    } else {
        plot <- plotTopCountDistributionByFactor(data, "MHC.allele", factor, top)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

## Plotting functions for CDR3 length
# For both CDR3a and CDR3b, x is the clone.id and y is the length of CDR3, colored by an inputted factor
plotCDR3SeqLength <- function(data, plotly = TRUE, filter.by = NULL, factor = "Epitope") {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    plot <- ggplot((data %>%
        pivot_longer(cols = c(CDR3a, CDR3b), names_to = "CDR3", values_to = "CDR3.seq")),
        aes(x = clone.id, y = str_length(CDR3.seq), color = !!sym(factor))) +
        geom_point() +
        facet_wrap(~CDR3, scales = "free_x") +
        theme_minimal() +
        theme(axis.text.x = element_blank(),
              strip.text = element_text(size = 8),
              legend.position = "none") +
        labs(title = "CDR3 Length",
                x = "Clone ID",
                y = "Length") +
        scale_color_viridis_d(option = "viridis")


    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}


## Plotting functions for full sequence clustering



## Plotting functions for CDR clustering



## Plotting functions for CDR3 clustering



## Plotting functions for full-receptor PCA plotting
# Common wrapper: Full seq, CDR3 seqs, V+J genes, MHC, color by epitope



## Generic functions
plotCountDistribution <- function(data, column) {
    # If top is not null, group all other values not in the top n counts into "Other" category
    data <- data %>%
        add_count(!!sym(column), name = "Count") %>%
        distinct(Epitope, .keep_all = TRUE)
    plot <- ggplot(data, aes(x = Count)) +
        geom_density() + 
        theme_minimal() +
        scale_x_log10() +
        labs(title = paste0(column, "Count Distribution"),
                x = column,
                y = "Density")
}


plotCountDistributionByFactor <- function(data, column, factor) {
    column.sym <- sym(column)
    column.n.sym <- sym(paste0(factor, " count"))
    factor.sym <- sym(factor)
    factor.n.sym <- sym(paste0(factor, " count"))
    data <- data %>% 
        select(!!column.sym, !!factor.sym) %>%
            add_count(!!column.sym, name = paste0(factor, " count")) %>%
            add_count(!!factor.sym, name = paste0(factor, " count")) %>%
            distinct(.keep_all = TRUE) %>%
            na.omit()
    plot <- ggplot(data, aes(x = !!column.n.sym, color = !!factor.sym)) +
        geom_density() + 
        theme_minimal() +
        scale_x_log10() +
        labs(title = paste0(column, "Count Distribution"),
                x = column,
                y = "Density") + 
        scale_color_viridis_d(option = "viridis")
}

plotTopCountDistributionByFactor <- function(data, column, factor, top = 10) {
    column.sym <- sym(column)
    column.n.sym <- sym(paste0(column, " count"))
    factor.sym <- sym(factor)
    factor.n.sym <- sym(paste0(factor, " count"))
    data <- data %>%
        select(!!column.sym, !!factor.sym) %>%
            add_count(!!column.sym, name = paste0(column, " count")) %>%
            add_count(!!factor.sym, name = paste0(factor, " count")) %>%
            distinct(.keep_all = TRUE) %>%
            na.omit()
    data <- data %>% 
        mutate(!!factor.sym := ifelse(!!factor.sym %in% 
        (data %>% 
            distinct(!!factor.sym, .keep_all = TRUE) %>% 
            slice_max(!!factor.n.sym, n = top) %>% 
            pull(!!factor.sym)), 
        !!factor.sym, "Other"))
    color.names <- data %>% distinct(!!factor.sym) %>% pull(!!factor.sym) 
    color.names <- color.names <- setdiff(color.names, "Other")
    color.factor <- c(viridis(length(color.names), option = "viridis"), "grey")
    names(color.factor) <- c(color.names, "Other")
    plot <- ggplot(data %>%
        mutate(!!factor.sym := factor(!!factor.sym, levels = c(color.names, "Other"))), aes(x = !!column.n.sym, color = !!factor.sym)) +
        geom_density() + 
        theme_minimal() +
        scale_x_log10() +
        labs(title = paste0(column, "Count Distribution"),
                x = column,
                y = "Density") + 
        scale_color_manual(values = color.factor)
}

plotComposition <- function(data, column) {
    data <- data %>%
            add_count(!!sym(column), name = "Count") 
    plot <- ggplot((data %>% distinct(!!sym(column), Count, .keep_all = TRUE)), 
        aes(x = !!sym(column), y = Count, fill = !!sym(column))) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "none") +
        labs(title = paste0(column, " Composition"),
                x = column,
                y = "Count") +
        scale_fill_viridis_d(option = "viridis")
}

plotCompositionByFactor <- function(data, column, factor) {
    column.sym <- sym(column)
    factor.sym <- sym(factor)
    data <- data %>%
            select(!!column.sym, !!factor.sym) %>%
            add_count(!!column.sym, name = "Count") %>%
            add_count(!!factor.sym, name = paste0(factor, " count")) %>%
            na.omit()
    plot <- ggplot(data, 
        aes(x = !!column.sym, y = `Count`, fill = !!factor.sym)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(vars(!!factor.sym), strip.position = "top", scales = "free_x", nrow = 1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "none") +
        labs(title = paste0(column, " Composition"),
                x = column,
                y = "Count") +
        scale_fill_viridis_d(option = "viridis")
}

plotTopCompositionByFactor <- function(data, column, factor, top = 10) {
    column.sym <- sym(column)
    factor.sym <- sym(factor)
    data <- data %>%
            select(!!column.sym, !!factor.sym) %>%
            add_count(!!column.sym, name = "Count") %>%
            add_count(!!factor.sym, name = paste0(factor, " count")) %>%
            distinct(.keep_all = TRUE) %>%
            na.omit()
    data <- data %>% 
        mutate(!!factor.sym := ifelse(!!factor.sym %in% 
        (data %>% 
            distinct(!!factor.sym, .keep_all = TRUE) %>% 
            slice_max(!!factor.sym, n = top) %>% 
            pull(!!factor.sym)), 
        !!factor.sym, "Other"))
    plot <- ggplot(data, 
        aes(x = !!column.sym, y = `Count`, fill = !!factor.sym)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(vars(!!factor.sym), strip.position = "top", scales = "free_x", nrow = 1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "none") +
        labs(title = paste0(column, " Composition"),
                x = column,
                y = "Count") +
        scale_fill_viridis_d(option = "viridis")
}

plotClusterResults <- function(data, column, chain, filter.by = NULL, color = NULL) {
    
    umap1.name <- paste0(column, "_UMAP_1")
    umap2.name <- paste0(column, "_UMAP_2")
    cluster.name <- paste0(column, "_cluster_0.5")
    if (is.null(color)) {
        color <- cluster.name
    }

    if (!umap1.name %in% colnames(data)) {
        data <- data %>%
            getCDRCluster(column, chain, filter.by = filter.by)
    } else if (umap1.name %in% colnames(data) & !is.null(filter.by)) {
        data <- data %>%
            filter(!!rlang::parse_expr(filter.by))
    }

    plot <- ggplot(data, aes_string(x = umap1.name, y = umap2.name, color = color)) +
        geom_point() +
        theme_minimal() +
        labs(title = paste0(column, " Clustering"),
                x = umap1.name,
                y = umap2.name, 
                color = paste0(column, " cluster")) +
        scale_color_viridis_d(option = "viridis")

    return(plot)

}

# Wrapper function for plotting with top approach
plotClusterResultsTop <- function(data, column, chain, filter.by = NULL, top = 20, top.col = NULL, color = NULL) {

    cluster.name <- paste0(column, "_cluster_0.5")

    if (is.null(color)) {
        color <- cluster.name
    }

    if (is.null(top.col)) {
        top.col <- color
    }

    data <- data %>%
        mutate(label = !!sym(color)) %>%
        add_count(label, name = "n")

    top_labels <- data %>%
        count(label) %>%
        top_n(top, n) %>%
        pull(label)

    data <- data %>%
        mutate(label = ifelse(label %in% top_labels, label, "Other"))

    labels <- unique(data$label)
    colors <- viridis(length(labels))
    names(colors) <- labels
    colors["Other"] <- "#717070"

    plot <- plotClusterResults(data, column, chain, filter.by, color = "label") +
        scale_color_manual(values = colors)

    return(plot)
}

# Wrapper function for plotting with threshold approach
plotClusterResultsThreshold <- function(data, column, chain, filter.by = NULL, threshold = 100, threshold.col = NULL, color = NULL) {
    cluster.name <- paste0(column, "_cluster_0.5")

    if (is.null(color)) {
        color <- cluster.name
    }

    if (is.null(threshold.col)) {
        threshold.col <- color
    }

    data <- data %>%
        rename(label = !!sym(color)) %>%
        add_count(label, name = "n") %>%
        mutate(label = ifelse(n >= threshold, label, "Other"))

    labels <- unique(data$label)
    colors <- viridis(length(labels), option ="mako")
    names(colors) <- labels
    colors["Other"] <- "#717070"

    plot <- plotClusterResults(data, column, chain, filter.by, color = "label") +
        scale_color_manual(values = colors)

    return(plot)
}

plotPCA <- function(data, columns, filter.by = NULL, color = "Epitope") {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    pca <- prcomp(data %>% select(all_of(columns)), scale = TRUE)
    pca.data <- data.frame(pca$x)
    pca.data <- cbind(data, pca.data)

    plot <- ggplot(pca.data, aes_string(x = "PC1", y = "PC2", color = color)) +
        geom_point() +
        theme_minimal() +
        labs(title = "PCA Plot",
                x = "PC1",
                y = "PC2",
                color = color) +
        scale_color_viridis_d(option = "viridis")

    return(plot)
}