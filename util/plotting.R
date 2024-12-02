

## Plotting functions for V/J genes

pivotVJData <- function(data) {

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
            filter(keep == `Gene Chain`)
    
    return(data)

}

plotVJAlleleCountDistribution <- function(data, plotly = TRUE, filter.by = NULL, group.by = NULL, top = NULL, option = "viridis", begin = 0, end = 1) {

    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    data <- data %>%
        pivotVJData() 

    if (is.null(top) && !is.null(group.by)) {
        plot <- plotCountDistributionByFactor(data, "Allele", group.by, option, begin, end) + 
            facet_wrap(~`Gene Segment`, scales = "free", ncol = 2)
    } else if (!is.null(top) && !is.null(group.by)) {
        plot <- plotTopCountDistributionByFactor(data, "Allele", group.by, top, option, begin, end) + 
            facet_wrap(~`Gene Segment`, scales = "free", ncol = 2)
    } else {
        plot <- plotCountDistribution(data, "Allele", option) + 
            facet_wrap(~`Gene Segment`, scales = "free", ncol = 2)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }

}

plotVJGeneCountDistribution <- function(data, plotly = TRUE, filter.by = NULL, group.by = NULL, top = NULL, option = "viridis", begin = 0, end = 1) {

    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    data <- data %>% 
        pivotVJData()

    if (is.null(top) && !is.null(group.by)) {
        plot <- plotCountDistributionByFactor(data, "Gene", group.by, option, begin, end) + 
            facet_wrap(~`Gene Segment`, scales = "free", ncol = 2)
    } else if (!is.null(top) && !is.null(group.by)) {
        plot <- plotTopCountDistributionByFactor(data, "Gene", group.by, top, option, begin, end) + 
            facet_wrap(~`Gene Segment`, scales = "free", ncol = 2)
    } else {
        plot <- plotCountDistribution(data, "Gene", option) + 
            facet_wrap(~`Gene Segment`, scales = "free", ncol = 2)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }

}

plotVJFamilyCountDistribution <- function(data, plotly = TRUE, filter.by = NULL, group.by = NULL, top = NULL, option = "viridis", begin = 0, end = 1) {

    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    data <- data %>% 
        pivotVJData()

    if (is.null(top) && !is.null(group.by)) {
        plot <- plotCountDistributionByFactor(data, "Gene Family", group.by, option, begin, end) + 
            facet_wrap(~`Gene Segment`, scales = "free", ncol = 2)
    } else if (!is.null(top) && !is.null(group.by)) {
        plot <- plotTopCountDistributionByFactor(data, "Gene Family", group.by, top, option, begin, end) + 
            facet_wrap(~`Gene Segment`, scales = "free", ncol = 2)
    } else {
        plot <- plotCountDistribution(data, "Gene Family", option) + 
            facet_wrap(~`Gene Segment`, scales = "free", ncol = 2)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }

}


plotVJAllele <- function(data, plotly = TRUE, filter.by = NULL, option = "viridis", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    data <- data %>%
        pivotVJData() %>%
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
        scale_fill_viridis_d(option = option, begin = begin, end = end)

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}


plotVJGene <- function(data, plotly = TRUE, filter.by = NULL, option = "viridis", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    data <- data %>%
            pivotVJData() %>%
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
        scale_fill_viridis_d(option = option, begin = begin, end = end)

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

plotVJFamily <- function(data, plotly = TRUE, filter.by = NULL, option = "viridis", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    data <- data %>% 
        pivotVJData()  %>% 
        add_count(`Gene Family`, name = "Count") %>%
        ungroup() %>% 
        mutate(`Gene Family` = as.character(`Gene Family`))

    plot <- ggplot(data, 
        aes(x = `Gene Family`, fill = `Gene Family`)) +
        geom_bar(position = "dodge") +
        facet_wrap(~`Gene Family`, strip.position = "bottom", scales = "free", nrow = 1) +
        facet_wrap(~`Gene Segment`, scales = "free") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              strip.text = element_text(size = 8),
              legend.position = "none") +
        labs(title = "V/J Gene Families",
                x = "Gene",
                y = "Count") +
        scale_fill_viridis_d(option = option, begin = begin, end = end)

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

# More generic VJ composition function
plotVJComposition <- function(data, plotly = TRUE, filter.by = NULL, group.by = NULL, option = "viridis", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    data <- data %>% 
        add_count(Allele, name = "Count") %>%
        distinct(.keep_all = TRUE)

    plot <- ggplot(data = data, 
        aes(x = Allele, y = Count, fill = Allele)) +
        geom_bar(stat = "identity") +
        facet_wrap(~`Gene Segment`, scales = "free") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              strip.text = element_text(size = 8),
              legend.position = "none") +
        labs(title = "V/J Alleles",
                x = "Allele",
                y = "Count") +
        scale_fill_viridis_d(option = option, begin = begin, end = end)

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}


# VJ sankey diagram
plotVJSankey <- function(data, filter.by = NULL, top = NULL, threshold = NULL) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    data <- data %>% 
        select(AV, AJ, BV, BJ)

    if (!is.null(top)) {
        top.av <- data %>% 
            count(AV) %>% 
            slice_max(order_by = n, n = top) %>%
            pull(AV)
        top.aj <- data %>% 
            count(AJ) %>%
            slice_max(order_by = n, n = top) %>%
            pull(AJ)
        top.bv <- data %>% 
            count(BV) %>%
            slice_max(order_by = n, n = top) %>%
            pull(BV)
        top.bj <- data %>% 
            count(BJ) %>%
            slice_max(order_by = n, n = top) %>%
            pull(BJ)
        data <- data %>% 
            mutate(AV = ifelse(AV %in% top.av, AV, "Other AV"),
                AJ = ifelse(AJ %in% top.aj, AJ, "Other AJ"),
                BV = ifelse(BV %in% top.bv, BV, "Other BV"),
                BJ = ifelse(BJ %in% top.bj, BJ, "Other AV"))
    } else if (!is.null(threshold)) {
        threshold.av <- data %>% 
            count(AV, name = "AV.count") %>%
            filter(AV.count >= threshold) %>%
            pull(AV)
        threshold.aj <- data %>%
            count(AJ, name = "AJ.count") %>%
            filter(AJ.count >= threshold) %>%
            pull(AJ)
        threshold.bv <- data %>%
            count(BV, name = "BV.count") %>%
            filter(BV.count >= threshold) %>%
            pull(BV)
        threshold.bj <- data %>%
            count(BJ, name = "BJ.count") %>%
            filter(BJ.count >= threshold) %>%
            pull(BJ)
        data <- data %>%
            mutate(AV = ifelse(AV %in% threshold.av, AV, "Other AV"),
                AJ = ifelse(AJ %in% threshold.aj, AJ, "Other AJ"),
                BV = ifelse(BV %in% threshold.bv, BV, "Other BV"),
                BJ = ifelse(BJ %in% threshold.bj, BJ, "Other BJ"))
    }

    vj.freq <- data %>% 
        count(AV, AJ, BV, BJ) %>%
        arrange(desc(n))

    nodes <- data.frame(name = unique(c(vj.freq$AV, vj.freq$AJ, vj.freq$BV, vj.freq$BJ)))
    links <- data.frame(source = match(vj.freq$AJ, nodes$name) - 1,
        target = match(vj.freq$AV, nodes$name) - 1,
        value = vj.freq$n,
        stringsAsFactors = FALSE)
    links <- rbind(links,
        data.frame(source = match(vj.freq$AV, nodes$name) - 1,
        target = match(vj.freq$BV, nodes$name) - 1,
        value = vj.freq$n,
        stringsAsFactors = FALSE))
    links <- rbind(links,
        data.frame(source = match(vj.freq$BV, nodes$name) - 1,
        target = match(vj.freq$BJ, nodes$name) - 1,
        value = vj.freq$n,
        stringsAsFactors = FALSE))
    nodes <- nodes %>%
        mutate(source = 0:(nrow(nodes) - 1))
    links <- links %>%
        left_join(nodes, join_by(source))
    links <- links %>%
        mutate(gene = name) %>%
        select(-name)

    sankeyplot <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source", 
        Target = "target", Value = "value", NodeID = "name", sinksRight = FALSE, LinkGroup = "gene")
}


# V and J gene heatmaps
# Plot V/J gene usage heatmaps
plotVJHeatmapNormalized <- function(data, gene.name) {
  allele.counts <- data %>%
    gather(key = "gene", value = "allele", AV, AJ, BV, BJ) %>%
    filter(gene == gene.name) %>%
    group_by(source, gene, allele) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(source) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    ungroup()

  heatmap.data <- dcast(allele.counts, source + gene ~ allele, value.var = "percentage", fill = 0)
  heatmap.matrix <- as.matrix(heatmap.data[,-c(1, 2)])
  rownames(heatmap.matrix) <- heatmap.data$source
  heatmap.data <- melt(heatmap.matrix)

  unique_alleles <- unique(heatmap.data$Var2)
#   ordered_alleles <- sortAlleles(unique_alleles)
  heatmap.data$Var2 <- factor(heatmap.data$Var2)

  heatmap <- ggplot(heatmap.data, aes(Var2, Var1, fill = value)) +
    geom_tile() +
     scale_fill_gradientn(colors = c("#d1e5f0", "#92c5de", "#f7fcb9", "#fee08b", "#fdae61", "#f46d43", "#d73027")) +
    labs(title = paste0("TR", gene.name, " Heatmap"),
         x = "Allele",
         y = "Sample",
         fill = "Normalized\npercentage") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

  return(heatmap)
}


# Wrapper for normalized heatmaps
makeVJHeatmapsNormalized <- function(data, column) {

  data <- data %>% 
    mutate(source = !!sym(column))
#   for (i in unique(data[[column]])) {
#     sample <- data[[column]] %>% mutate(source = data.names[i])
#     data.full <- rbind(data.full, sample)
#   }

    AV <- plotVJHeatmapNormalized(data, "AV")
    AJ <- plotVJHeatmapNormalized(data, "AJ")
    BV <- plotVJHeatmapNormalized(data, "BV")
    BJ <- plotVJHeatmapNormalized(data, "BJ")

    return(list(AV, AJ, BV, BJ))
}



## Plotting functions for epitope
plotEpitopeCountDistribution <- function(data, plotly = TRUE, filter.by = NULL, group.by = NULL, top = NULL, option = "turbo", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    if (is.null(top) && !is.null(group.by)) {
        plot <- plotCountDistributionByFactor(data, "Epitope", group.by, option, begin, end)
    } else if (!is.null(top) && !is.null(group.by)) {
        plot <- plotTopCountDistributionByFactor(data, "Epitope", group.by, top, option, begin, end)
    } else {
        plot <- plotCountDistribution(data, "Epitope", option)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}


plotEpitopeComposition <- function(data, plotly = TRUE, filter.by = NULL, group.by = NULL, top = NULL, option = "turbo", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    if (is.null(top) && !is.null(group.by)) {
        plot <- plotCompositionByFactor(data, "Epitope", group.by, option, begin, end)
    } else if (!is.null(top) && !is.null(group.by)) {
        plot <- plotTopCompositionByFactor(data, "Epitope", group.by, top, option, begin, end)
    } else {
        plot <- plotComposition(data, "Epitope", option, begin, end)
    }

    if (plotly) {
        if (is.null(filter.by) && is.null(group.by)) {
            plot <- plot + theme(axis.text.x = element_blank())
        } 
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}


plotEpitopeSpeciesComposition <- function(data, plotly = TRUE, filter.by = NULL, group.by = NULL, top = NULL, option = "turbo", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    if (is.null(top) && !is.null(group.by)) {
        plot <- plotCompositionByFactor(data, "Epitope.species", group.by, option, begin, end)
    } else if (!is.null(top) && !is.null(group.by)) {
        plot <- plotTopCompositionFactor(data, "Epitope.species", group.by, top, option, begin = begin, end = end)
    } else {
        plot <- plotComposition(data, "Epitope.species", option, begin = begin, end = end)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

## Plotting functions for MHC
plotMHCAlleleCountDistribution <- function(data, plotly = TRUE, filter.by = NULL, group.by = NULL, top = NULL, option = "plasma", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    if (is.null(top) && !is.null(group.by)) {
        plot <- plotCountDistributionByFactor(data, "MHC.allele", group.by, option, begin = begin, end = end) +
            scale_fill_viridis_d(option = option, begin = 1, end = 0.1)
    } else if (!is.null(top) && !is.null(group.by)) {
        plot <- plotTopCountDistributionByFactor(data, "MHC.allele", group.by, top, option, begin = begin, end = end) +
            scale_fill_viridis_d(option = option, begin = 1, end = 0.1)
    } else {
        plot <- plotCountDistribution(data, "MHC.allele", option) 
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

plotMHCLocusCountDistribution <- function(data, plotly = TRUE, filter.by = NULL, group.by = NULL, top = NULL, option = "plasma", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    if (is.null(top) && !is.null(group.by)) {
        plot <- plotCountDistributionByFactor(data, "MHC.locus", group.by, option, begin = begin, end = end)
    } else if (!is.null(top) && !is.null(group.by)) {
        plot <- plotTopCountDistributionByFactor(data, "MHC.locus", group.by, top, option, begin = begin, end = end)
    } else {
        plot <- plotCountDistribution(data, "MHC.locus", option)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}


plotMHCCountDistribution <- function(data, plotly = TRUE, filter.by = NULL, group.by = NULL, top = NULL, option = "plasma", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    if (is.null(top) && !is.null(group.by)) {
        plot <- plotCountDistributionByFactor(data, "MHC", group.by, option, begin = begin, end = end)
    } else if (!is.null(top) && !is.null(group.by)) {
        plot <- plotTopCountDistributionByFactor(data, "MHC", group.by, top, option, begin = begin, end = end)
    } else {
        plot <- plotCountDistribution(data, "MHC", option)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

plotMHCComposition <- function(data, plotly = TRUE, filter.by = NULL, group.by = NULL, top = NULL, option = "plasma", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }
    
    if (is.null(top) && !is.null(group.by)) {
        plot <- plotCompositionByFactor(data, "MHC", group.by, option, begin = begin, end = end)
    } else if (!is.null(top) && !is.null(group.by)) {
        plot <- plotTopCompositionByFactor(data, "MHC", group.by, top, option, begin = begin, end = end)
    } else {
        plot <- plotComposition(data, "MHC", option, begin = begin, end = end)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}


plotMHCLocusComposition <- function(data, plotly = TRUE, filter.by = NULL, group.by = NULL, top = NULL, option = "plasma", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    if (is.null(top) && !is.null(group.by)) {
        plot <- plotCompositionByFactor(data, "MHC.locus", group.by, option, begin = begin, end = end)
    } else if (!is.null(top) && !is.null(group.by)) {
        plot <- plotTopCompositionByFactor(data, "MHC.locus", group.by, top, option, begin = begin, end = end)
    } else {
        plot <- plotComposition(data, "MHC.locus", option, begin = begin, end = end)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}


plotMHCAlleleComposition <- function(data, plotly = TRUE, filter.by = NULL, group.by = NULL, top = NULL, option = "plasma", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    if (is.null(top) && !is.null(group.by)) {
        plot <- plotCompositionByFactor(data, "MHC.allele", group.by, option, begin = begin, end = end)
    } else if (!is.null(top) && !is.null(group.by)) {
        plot <- plotTopCompositionByFactor(data, "MHC.allele", group.by, top, option, begin = begin, end = end)
    } else {
        plot <- plotComposition(data, "MHC.allele", option, begin = begin, end = end)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}


## Plotting functions for CDR3 length
# For both CDR3a and CDR3b, x is the clone.id and y is the length of CDR3, colored by an inputted group.by
plotCDR3SeqLength <- function(data, color.by = "Epitope", plotly = TRUE, filter.by = NULL, top = NULL, option = "viridis", begin = 0, end = 1) {

    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

  
    plot <- ggplot((data %>%
        pivot_longer(cols = c(CDR3a, CDR3b), names_to = "CDR3.chain", values_to = "CDR3.seq")),
        aes(x = clone.id, y = str_length(CDR3.seq), color = !!sym(color.by))) +
        geom_point() +
        facet_wrap(~CDR3.chain, scales = "free_x") +
        theme_minimal() +
        theme(axis.text.x = element_blank(),
              strip.text = element_text(size = 8),
              legend.position = "bottom") +
        labs(title = "CDR3 Length",
                x = "Clone ID",
                y = "Length") +
        scale_color_viridis_d(option = option, begin = begin, end = end)

    if (!is.null(top)) {
        data <- data %>%
            mutate(label = !!sym(color.by)) %>%
            add_count(label, name = "n")

        top.labels <- data %>%
            count(label) %>%
            slice_max(n = top, order_by = n) %>%
            pull(label)

        data <- data %>%
            mutate(label = ifelse(label %in% top.labels, label, "Other"))

        labels <- unique(data$label)
        colors <- viridis(length(labels), option = option, begin = begin, end = end)
        names(colors) <- labels
        colors["Other"] <- "#717070"
        plot <- plot + scale_color_manual(values = colors)
    }


    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

# Density plot of max similarity in CDR3a and CDR3b
plotCDR3SimilarityDistribution <- function(data, plotly = TRUE, filter.by = NULL, top = NULL, option = "mako", begin = 0.2, end = 0.8) {

    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    plot <- ggplot(data, aes(x = CDR3a.max.similarity, fill = "CDR3a")) +
        geom_density(alpha = 0.5) +
        geom_density(aes(x = CDR3b.max.similarity, fill = "CDR3b"), alpha = 0.5) +
        scale_fill_viridis_d(option = option, begin = begin, end = end) +
        theme_minimal() +
        labs(title = "Max CDR3 Similarity Density",
                x = "Max Similarity",
                y = "Density")

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}

# Box plot of max similarity in CDR3a and CDR3b
plotCDR3SimilarityByFactor <- function(data, plotly = TRUE, filter.by = NULL, top = NULL, group.by = NULL, option = "mako", begin = 0.2, end = 0.8) {

    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    group.by.sym <- sym(group.by)

    data <- data %>% mutate(!!group.by.sym := ifelse(!!group.by.sym %in% (data %>% count(!!group.by.sym) %>% 
        slice_max(order_by = n, n = 10) %>% pull(!!group.by.sym )), !!group.by.sym, "Other"))

    data <- data %>%
        pivot_longer(cols = c(CDR3a.max.similarity, CDR3b.max.similarity), names_to = "CDR3.chain", values_to = "CDR3.similarity")

    plot <- ggplot(data, aes_string(x = group.by, y = "CDR3.similarity", fill = "CDR3.chain")) +
        geom_boxplot() +
        scale_fill_viridis_d(option = option, begin=begin, end=end, alpha=0.5) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title.x = element_blank(),
              legend.position = "none") +
        labs(title = "CDR3a Similarity",
                x = group.by,
                y = "Max. similarity")

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}


## Plotting functions for full sequence clustering
plotFullSeqClusterResults <- function(data, plotly = TRUE, filter.by = NULL, color = NULL, top = NULL, top.col = NULL, threshold = NULL, threshold.col = NULL, highlight = NULL, option = "viridis", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    if (!is.null(top)) {
        plot <- plotClusterResultsTop(data, "full.seq", "AB", top = top, top.col = top.col, option = option, begin = begin, end = end)
    } else if (!is.null(threshold)) {
        plot <- plotClusterResultsThreshold(data, "full.seq", "AB", threshold = threshold, threshold.col = threshold.col, option = option, begin = begin, end = end)
    } else if (!is.null(highlight)) {
        plot <- plotClusterResultsHighlight(data, "full.seq", "AB", filter.by = filter.by, color = color, highlight = highlight, option = option, begin = begin, end = end)
    } else {
        plot <- plotClusterResults(data, "full.seq", "AB", filter.by = filter.by, color = color, option = option, begin = begin, end = end)
    }

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}


## Plotting functions for CDR clustering
plotCDRSeqClusterResults <- function(data, combine.chains = TRUE, plotly = TRUE, filter.by = NULL, color = NULL, top = NULL, top.col = NULL, threshold = NULL, threshold.col = NULL, highlight = NULL, option = "viridis", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    if (combine.chains) {
        data <- data %>%
            mutate(CDR.seq = paste0(CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b))

        if (!is.null(top)) {
            plot <- plotClusterResultsTop(data, "CDR.seq", "AB", top = top, top.col = top.col, option = option, begin = begin, end = end)
        } else if (!is.null(threshold)) {
            plot <- plotClusterResultsThreshold(data, "CDR.seq", "AB", threshold = threshold, threshold.col = threshold.col, option = option, begin = begin, end = end)
        } else if (!is.null(highlight)) {
            plot <- plotClusterResultsHighlight(data, "CDR.seq", "AB", filter.by = filter.by, color = color, highlight = highlight, option = option, begin = begin, end = end)
        } else {
            plot <- plotClusterResults(data, "CDR.seq", "AB", filter.by = filter.by, color = color, option = option, begin = begin, end = end)
        }

        if (plotly) {
            plot <- ggplotly(plot)
            return(plot)
        } else {
            return(plot)
        }

    } else {
        data <- data %>% 
            mutate(CDRa.seq = paste0(CDR1a, CDR2a, CDR2.5a, CDR3a),
                   CDRb.seq = paste0(CDR1b, CDR2b, CDR2.5b, CDR3b))

        if (!is.null(top)) {
            plota <- plotClusterResultsTop(data, "CDRa.seq", "A", top = top, top.col = top.col, option = option, begin = begin, end = end)
            plotb <- plotClusterResultsTop(data, "CDRb.seq", "B", top = top, top.col = top.col, option = option, begin = begin, end = end) + 
                labs(title = "CDR3 Clustering")
        } else if (!is.null(threshold)) {
            plota <- plotClusterResultsThreshold(data, "CDRa.seq", "A", threshold = threshold, threshold.col = threshold.col, option = option, begin = begin, end = end)
            plotb <- plotClusterResultsThreshold(data, "CDRb.seq", "B", threshold = threshold, threshold.col = threshold.col, option = option, begin = begin, end = end) + 
                labs(title = "CDR3 Clustering")
        } else if (!is.null(highlight)) {
            plota <- plotClusterResultsHighlight(data, "CDRa.seq", "A", filter.by = filter.by, color = color, highlight = highlight, option = option, begin = begin, end = end)
            plotb <- plotClusterResultsHighlight(data, "CDRb.seq", "B", filter.by = filter.by, color = color, highlight = highlight, option = option, begin = begin, end = end) + 
                labs(title = "CDR3 Clustering")
        } else {
            plota <- plotClusterResults(data, "CDRa.seq", "A", filter.by = filter.by, color = color, option = option, begin = begin, end = end)
            plotb <- plotClusterResults(data, "CDRb.seq", "B", filter.by = filter.by, color = color, option = option, begin = begin, end = end) + 
                labs(title = "CDR3 Clustering")
        }

        if (plotly) {
            plot <- subplot(plota, plotb, nrows = 1)
            return(plot)
        } else {
            plot <- ggarrange(plota, plotb, nrow = 1, common.legend = TRUE)
            return(plot)
        }
    }

}


## Plotting functions for CDR3 clustering
plotCDR3SeqClusterResults <- function(data, combine.chains = TRUE, plotly = TRUE, filter.by = NULL, color = NULL, top = NULL, top.col = NULL, threshold = NULL, threshold.col = NULL, highlight = NULL, option = "viridis", begin = 0, end = 1) {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    if (combine.chains) {
        data <- data %>%
            mutate(CDR3 = paste0(CDR3a, CDR3b))

        if (!is.null(top)) {
            plot <- plotClusterResultsTop(data, "CDR3", "AB", top = top, top.col = top.col, option = option, begin = begin, end = end)
        } else if (!is.null(threshold)) {
            plot <- plotClusterResultsThreshold(data, "CDR3", "AB", threshold = threshold, threshold.col = threshold.col, option = option, begin = begin, end = end)
        } else if (!is.null(highlight)) {
            plot <- plotClusterResultsHighlight(data, "CDR3", "AB", filter.by = filter.by, color = color, highlight = highlight, option = option, begin = begin, end = end)
        } else {
            plot <- plotClusterResults(data, "CDR3", "AB", filter.by = filter.by, color = color, option = option, begin = begin, end = end)
        }

        if (plotly) {
            plot <- ggplotly(plot)
            return(plot)
        } else {
            return(plot)
        }

    } else {

        if (!is.null(top)) {
            plota <- plotClusterResultsTop(data, "CDR3a", "A", top = top, top.col = top.col, option = option, begin = begin, end = end)
            plotb <- plotClusterResultsTop(data, "CDR3b", "B", top = top, top.col = top.col, option = option, begin = begin, end = end) + 
                labs(title = "CDR3 Clustering")
        } else if (!is.null(threshold)) {
            plota <- plotClusterResultsThreshold(data, "CDR3a", "A", threshold = threshold, threshold.col = threshold.col, option = option, begin = begin, end = end)
            plotb <- plotClusterResultsThreshold(data, "CDR3b", "B", threshold = threshold, threshold.col = threshold.col, option = option, begin = begin, end = end) + 
                labs(title = "CDR3 Clustering")
        } else if (!is.null(highlight)) {
            plota <- plotClusterResultsHighlight(data, "CDR3a", "A", filter.by = filter.by, color = color, highlight = highlight, option = option, begin = begin, end = end)
            plotb <- plotClusterResultsHighlight(data, "CDR3b", "B", filter.by = filter.by, color = color, highlight = highlight, option = option, begin = begin, end = end) + 
                labs(title = "CDR3 Clustering")
        } else {
            plota <- plotClusterResults(data, "CDR3a", "A", filter.by = filter.by, color = color, option = option, begin = begin, end = end)
            plotb <- plotClusterResults(data, "CDR3b", "B", filter.by = filter.by, color = color, option = option, begin = begin, end = end) + 
                labs(title = "CDR3 Clustering")
        }

        if (plotly) {
            plot <- subplot(plota, plotb, nrows = 1)
            return(plot)
        } else {
            plot <- ggarrange(plota, plotb, nrow = 1, common.legend = TRUE)
            return(plot)
        }
    }
}


## Plotting functions for full-receptor PCA plotting
# Common wrapper: Full seq, CDR seqs, V+J genes, MHC, color by epitope (given PCA deals with correlated variables, it's ok to include full seq and all seq components as well as VJ gene + gene family data)
plotReceptorPCA <- function(data, plotly = TRUE, filter.by = NULL, color = "Epitope") {
    if (!is.null(filter.by)) {
        filter.expression <- rlang::parse_expr(filter.by)
        data <- data %>% filter(!!filter.expression)
    }

    columns <- c("full.seq", "CDR1a", "CDR2a", "CDR2.5a", "CDR3a", "CDR1b", "CDR2b", "CDR2.5b", "CDR3b", 
        "AV", "AJ", "BV", "BJ", "AV.family", "AJ.family", "BV.family", "BJ.family",
        "MHC", "MHC.locus")

    plot <- plotPCA(data, columns = columns, color = color)

    if (plotly) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}


## Generic functions
plotCountDistribution <- function(data, column, option = "viridis") {
    # If top is not null, group all other values not in the top n counts into "Other" category
    data <- data %>%
        add_count(!!sym(column), name = "Count") %>%
        distinct(!!sym(column), .keep_all = TRUE)
    plot <- ggplot(data, aes(x = Count)) +
        geom_density(color = viridis(1, option = option, begin = begin, end = end), fill = viridis(1, option = option, begin = begin, end = end), alpha = 0.7) + 
        theme_minimal() +
        scale_x_log10() +
        labs(title = paste0(column, " Count Distribution"),
                x = "Count",
                y = "Density")
}


plotCountDistributionByFactor <- function(data, column, group.by, option = "viridis", begin = 0, end = 1) {
    column.sym <- sym(column)
    column.n.sym <- sym(paste0(group.by, " count"))
    group.by.sym <- sym(group.by)
    group.by.n.sym <- sym(paste0(group.by, " count"))
    data <- data %>% 
        select(!!column.sym, !!group.by.sym) %>%
            add_count(!!column.sym, name = paste0(group.by, " count")) %>%
            add_count(!!group.by.sym, name = paste0(group.by, " count")) %>%
            distinct(!!column.sym, !!group.by.sym, .keep_all = TRUE) %>%
            na.omit()
    plot <- ggplot(data, aes(x = !!column.n.sym, color = !!group.by.sym, fill = !!group.by.sym)) +
        geom_density() + 
        theme_minimal() +
        scale_x_log10() +
        labs(title = paste0(column, " Count Distribution"),
                x = "Count",
                y = "Density") + 
        scale_color_viridis_d(option = option)
}

plotTopCountDistributionByFactor <- function(data, column, group.by, top = 10, option = "viridis", begin = 0, end = 1) {
    column.sym <- sym(column)
    column.n.sym <- sym(paste0(column, " count"))
    group.by.sym <- sym(group.by)
    group.by.n.sym <- sym(paste0(group.by, " count"))
    data <- data %>%
        select(!!column.sym, !!group.by.sym) %>%
            add_count(!!column.sym, name = paste0(column, " count")) %>%
            add_count(!!group.by.sym, name = paste0(group.by, " count")) %>%
            distinct(!!column.sym, !!group.by.sym, .keep_all = TRUE) %>%
            na.omit()
    data <- data %>% 
        mutate(!!group.by.sym := ifelse(!!group.by.sym %in% 
        (data %>% 
            distinct(!!group.by.sym, .keep_all = TRUE) %>% 
            slice_max(order_by = !!group.by.n.sym, n = top) %>% 
            pull(!!group.by.sym)), 
        !!group.by.sym, "Other"))
    color.names <- data %>% distinct(!!group.by.sym) %>% pull(!!group.by.sym) 
    color.names <- color.names <- setdiff(color.names, "Other")
    color.group.by <- c(viridis(length(color.names), option = option, begin = begin, end = end), "grey")
    names(color.group.by) <- c(color.names, "Other")
    plot <- ggplot(data %>%
        mutate(!!group.by.sym := factor(!!group.by.sym, levels = c(color.names, "Other"))), aes(x = !!column.n.sym, color = !!group.by.sym)) +
        geom_density() + 
        theme_minimal() +
        scale_x_log10() +
        labs(title = paste0(column, " Count Distribution"),
                x = "Count",
                y = "Density") + 
        scale_color_manual(values = color.group.by)
}

plotComposition <- function(data, column, option = "viridis", begin = 0, end = 1) {
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
        scale_fill_viridis_d(option = option, begin = begin, end = end)
}

plotCompositionByFactor <- function(data, column, group.by, option = "viridis", begin = 0, end = 1) {
    column.sym <- sym(column)
    group.by.sym <- sym(group.by)
    data <- data %>%
            select(!!column.sym, !!group.by.sym) %>%
            add_count(!!column.sym, name = "Count") %>%
            add_count(!!group.by.sym, name = paste0(group.by, " count")) %>%
            na.omit()
    plot <- ggplot(data, 
        aes(x = !!column.sym, y = `Count`, fill = !!group.by.sym)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(vars(!!group.by.sym), strip.position = "top", scales = "free_x", nrow = 1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "none") +
        labs(title = paste0(column, " Composition"),
                x = column,
                y = "Count") +
        scale_fill_viridis_d(option = option, begin = begin, end = end)
}

plotTopCompositionByFactor <- function(data, column, group.by, top = 10, option = "viridis", begin = 0, end = 1) {
    column.sym <- sym(column)
    column.n.sym <- sym(paste0(column, " count"))
    group.by.sym <- sym(group.by)
    group.by.n.sym <- sym(paste0(group.by, " count"))
    data <- data %>%
            select(!!column.sym, !!group.by.sym) %>%
            add_count(!!column.sym, name = "Count") %>%
            add_count(!!group.by.sym, name = paste0(group.by, " count")) %>%
            distinct(!!column.sym, !!group.by.sym, .keep_all = TRUE) %>%
            na.omit()
    data <- data %>% 
        filter(!!group.by.sym %in% (data %>% 
            distinct(!!group.by.sym, .keep_all = TRUE) %>% 
            slice_max(order_by = !!group.by.n.sym, n = top) %>% 
            pull(!!group.by.sym)))
    plot <- ggplot(data, 
        aes(x = !!column.sym, y = `Count`, fill = !!group.by.sym)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(vars(!!group.by.sym), strip.position = "top", scales = "free_x", nrow = 1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "none") +
        labs(title = paste0(column, " Composition"),
                x = column,
                y = "Count") +
        scale_fill_viridis_d(option = option, begin = begin, end = end)
}

plotClusterResults <- function(data, column, chain, filter.by = NULL, color = NULL, option = "viridis", begin = 0, end = 1) {
    
    umap1.name <- paste0(column, "_UMAP_1")
    umap2.name <- paste0(column, "_UMAP_2")
    cluster.name <- paste0(column, "_cluster_0.5")
    if (is.null(color)) {
        color <- cluster.name
    }
    data <- data %>% mutate(!!sym(color) := factor(!!sym(color)))


    if (!umap1.name %in% colnames(data)) {
        data <- data %>%
            getSeqCluster(column, chain, filter.by = filter.by)
    } else if (umap1.name %in% colnames(data) && !is.null(filter.by)) {
        data <- data %>%
            select(-c(umap1.name, umap2.name, cluster.name)) %>%
            getSeqCluster(column, chain, filter.by = filter.by)
    }

    plot <- ggplot(data, aes_string(x = umap1.name, y = umap2.name, color = color)) +
        geom_point(alpha = 0.7) +
        theme_minimal() +
        theme(legend.position = "bottom") +
        labs(title = paste0(column, " Clustering"),
                x = umap1.name,
                y = umap2.name, 
                color = paste0(column, " cluster")) +
        scale_color_viridis_d(option = option, begin = begin, end = end)  

    if (color != cluster.name) {
        plot <- plot + labs(color = color)
    }

    return(plot)

}

# Wrapper function for plotting with top approach
plotClusterResultsTop <- function(data, column, chain, filter.by = NULL, top = 20, top.col = NULL, option = "viridis", begin = 0, end = 1) {

    cluster.name <- paste0(column, "_cluster_0.5")

    if (is.null(top.col)) {
        top.col <- cluster.name
        data <- data %>% mutate(!!sym(color) := factor(!!sym(color)))
    } 

    if (!is.null(filter.by)) {
        data <- data %>%
            filter(!!rlang::parse_expr(filter.by))
    }


    data <- data %>%
        mutate(label = !!sym(top.col)) %>%
        add_count(label, name = "n")

    top.labels <- data %>%
        count(label) %>%
        slice_max(n = top, order_by = n) %>%
        pull(label)

    data <- data %>%
        mutate(label = ifelse(label %in% top.labels, label, "Other"))

    labels <- unique(data$label)
    colors <- viridis(length(labels), option = option, begin = begin, end = end)
    names(colors) <- labels
    colors["Other"] <- "#717070"

    plot <- plotClusterResults(data, column, chain, filter.by, color = "label") +
        scale_color_manual(values = colors) + 
        labs(color = top.col)

    return(plot)
}

# Wrapper function for plotting with threshold approach
plotClusterResultsThreshold <- function(data, column, chain, filter.by = NULL, threshold = 100, threshold.col = NULL, option = option, begin = 0, end = 1) {
    cluster.name <- paste0(column, "_cluster_0.5")

    if (is.null(threshold.col)) {
        threshold.col <- cluster.name
        data <- data %>% mutate(!!sym(color) := factor(!!sym(color)))
    }   

    if (!is.null(filter.by)) {
        data <- data %>%
            filter(!!rlang::parse_expr(filter.by))
    }

    data <- data %>%
        mutate(label = !!sym(threshold.col)) %>%
        add_count(label, name = "n") %>%
        mutate(label = ifelse(n >= threshold, label, "Other"))

    labels <- unique(data$label)
    colors <- viridis(length(labels), option = option, begin = begin, end = end)
    names(colors) <- labels
    colors["Other"] <- "#717070"

    plot <- plotClusterResults(data, column, chain, filter.by, color = "label") +
        scale_color_manual(values = colors) + 
        labs(color = threshold.col)

    return(plot)
}

# Wrapper function for plotting with specific values overlaid on all data (different from filtering as it shows all values in original clustering in grey, but uses the same expression syntax)
plotClusterResultsHighlight <- function(data, column, chain, filter.by = NULL, highlight, color, option = "viridis", begin = 0, end = 1) {

    if (!is.null(filter.by)) {
        data <- data %>%
            filter(!!rlang::parse_expr(filter.by))
    }

    highlight.labels <- data %>% 
        filter(!!rlang::parse_expr(highlight)) %>%
        mutate(label = !!sym(color)) %>%
        pull(!!sym(color))

    data <- data %>% 
        mutate(label = !!sym(color)) %>%
        mutate(label = ifelse(label %in% highlight.labels, label, "Other"))

    highlight <- unique(highlight.labels)
    colors <- viridis(length(highlight), option = option, begin = begin, end = end)
    names(colors) <- highlight
    colors["Other"] <- "#717070"

    plot <- plotClusterResults(data, column, chain, filter.by, color = "label") +
        scale_color_manual(values = colors) + 
        labs(color = color)

    return(plot)

}

plotPCA <- function(data, columns, filter.by = NULL, color = "Epitope", option = "viridis") {
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
        scale_color_viridis_d(option = option)

    return(plot)
}