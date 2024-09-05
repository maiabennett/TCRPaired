#########################
# Visualization methods #
#########################

# Plot clonotype distributions
plotClonotypeDistributions <- function(data, data.name) {
  data <- data %>% 
    mutate(clonotype = as.numeric(gsub("clonotype", "", clonotype)),
      clonotype.order = as.numeric(factor(clonotype)))
  data %>%
    ggplot(aes(x = clonotype.order, fill = as.factor(clonotype.order))) +
      geom_bar() +
      labs(title = paste0(data.name, " Clonotype Distribution"),
        x = "Clonotype",
        y = "Frequency") +
      scale_fill_viridis(option = "F", discrete = TRUE, begin = 0.8, end = 0.2) +
      theme(legend.position = "none")
}

# Plot epitope distributions
plotEpitopeDistributions <- function(data, fill, data.name) {
  data <- data %>% 
    mutate(Epitope = as.factor(Epitope))
  data %>%
    ggplot(aes(x = Epitope, fill = fill)) +
      geom_bar() +
      labs(title = paste0(data.name, " Epitope Distribution"),
        x = "Epitope",
        y = "Frequency") +
      scale_fill_viridis(option = "F", discrete = TRUE, begin = 0.8, end = 0.2) +
      theme(legend.position = "none")
}

# Plot epitope distributions as a stacked bar chart per some factor
plotEpitopeFactorDist <- function(data, x, colors) {
  data[[x]] <- factor(data[[x]], levels = unique(data[[x]]))
  data %>%
    #na.omit() %>%
    ggplot(aes(x = .data[[x]], fill = Epitope)) +
    geom_bar(position = "stack") +
    labs(title = paste0("Predicted epitopes per ", x),
      x = x,
      y = "Frequency") +
    scale_fill_manual(values = colors) +
    theme(legend.position = "bottom")
}

# Plot sequence similarity between receptors, colored by epitope prediction
plotSeqSimilarity <- function(data, data.name, column = "sTCR.epitope", colors) {
  data[[column]] <- factor(data[[column]], levels = unique(data[[column]]))
  data %>%
    ggplot(aes_string(x = "CDR3a.similarity", y = "CDR3b.similarity", color = column)) +
    geom_point(aes(size = 4)) +
    labs(title = paste0("Similarity between experimental and reference CDR3 sequences for ", data.name),
      x = "CDR3a Similarity",
      y = "CDR3b Similarity",
      color = "Predicted epitope") +
    scale_color_manual(values = colors) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme(legend.position = "bottom")
}

# Plot MHC group distributions


# Plot V/J gene family distributions


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
  ordered_alleles <- sortAlleles(unique_alleles)
  heatmap.data$Var2 <- factor(heatmap.data$Var2, levels = ordered_alleles)

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

# Gene co-usage for a single sample
plotGeneVSGeneHeatmap <- function(data, gene1, gene2, sample.name) {
  gene.counts <- data %>%
    count(!!sym(gene1), !!sym(gene2)) %>%
    complete(!!sym(gene1), !!sym(gene2), fill = list(n = 0))

  unique_alleles <- unique(gene.counts[[gene1]])
  ordered_alleles <- sortAlleles(unique_alleles)
  gene.counts[[gene1]] <- factor(gene.counts[[gene1]], levels = ordered_alleles)

  heatmap <- ggplot(gene.counts, aes_string(gene1, gene2, fill = "n")) +
    geom_tile() +
    scale_fill_gradientn(colors = c("#d1e5f0", "#92c5de", "#f7fcb9", "#fee08b", "#fdae61", "#f46d43", "#d73027"),
      limits = c(0, max(gene.counts$n, na.rm = TRUE))) +
    labs(title = paste0(gene1, " vs ", gene2, " Heatmap, ", sample.name),
         x = gene1,
         y = gene2, 
         fill = "Count") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  return(heatmap)
}

# Wrapper for normalized heatmaps, for comparing V/J gene usage between two samples (will need expanding as samples added)
makeVJHeatmapsNormalized <- function(data, data.names) {
  #data1 <- data1 %>% mutate(source = data.name1)
  #data2 <- data2 %>% mutate(source = data.name2)
  data.full <- data.frame()
  for (i in seq_along(data.names)) {
    sample <- data[[i]] %>% mutate(source = data.names[i])
    data.full <- rbind(data.full, sample)
  }

  AV <- plotVJHeatmapNormalized(data.full, "AV")
  AJ <- plotVJHeatmapNormalized(data.full, "AJ")
  BV <- plotVJHeatmapNormalized(data.full, "BV")
  BJ <- plotVJHeatmapNormalized(data.full, "BJ")

  return(list(AV, AJ, BV, BJ))
}



## General functions to assist in plotting
# Sorting V and J alleles for plots 
sortAlleles <- function(alleles) {
  alleles <- as.character(alleles)
  
  numeric_parts <- sapply(alleles, function(x) {
    parts <- gsub("[A-Z]", "", unlist(strsplit(x, "-")))
    parts <- unlist(strsplit(parts, "/"))
    as.numeric(parts[1])
  })
  
  order_indices <- order(numeric_parts)
  alleles[order_indices]
}

# Counting the number of times a labeling (overlay) value appears in the data to facilitate top n + "other" labeling
countLabelValues <- function(data, col) {
  count.label <- data %>%
    count(!!sym(col)) %>%
    rename_with(~"label", .cols = all_of(col))
  return(count.label)
}

# Relevel the data during plotting to whichever levels are specified
# Common use is to order clonotypes and clusters numerically, or to relevel T cell types
relevelPlot <- function(data, count.label, levels) {
  levels_mapping <- setNames(nm = levels, levels)
  for (level in levels) {
    if (level %in% count.label$label) {
      levels_mapping[level] <- paste0(level, "\nn = ", (count.label %>% filter(label == level) %>% select(n)))
    } 
  }
  levels <- levels_mapping
  data$label <- factor(data$label, levels = levels)
  return(data)
}


## Method set for plotting clusters with optional overlays ##

# Makes individual CDR cluster plots, with magnitude-aware point sizes for data where one point (identical CDR3) may represent multiple receptors
plotCDRCluster <- function(data, count.label, col, top, point.size, umap_x, umap_y, title, levels = NULL, method = FALSE) {

  if (!any(data$label == "Did not meet threshold", na.rm = TRUE)) {
    data <- data %>%
      mutate(label = ifelse(label != "Other", paste0(label, "\nn = ", n), "Other"))
  }

  else  {
    data <- data %>% 
      mutate(label = ifelse(label != "Did not meet threshold", label, "Did not meet threshold"),
        n = ifelse(label != "Did not meet threshold", n, 1))
  }

  if (!is.null(levels)) {
    data <- relevelPlot(data, count.label, levels)
  } else {
    data$label <- factor(data$label, levels = c(setdiff(unique(data$label), c("Other", NA)), "Other", NA))
  }

  if (!any(names(data) == "color")) {
    data$color <- viridis(length(unique(data$label)))
  } 

  color.mapping <- setNames(data$color, data$label)

  if (method) {
    clustplot <- data %>%
      filter(!label %in% c("Other", "Did not meet threshold") & !is.na(label)) %>% 
    ggplot(aes_string(x = umap_x, y = umap_y, color = "label")) +
      geom_point(data = data %>% filter(label %in% c("Other", "Did not meet threshold") | is.na(label)), 
        aes(x = !!sym(umap_x), y = !!sym(umap_y), color = label, size = 0.25),
        alpha = 0.4) +
      geom_point(aes(size = n, shape = Method), 
        alpha = 0.8) +
      labs(title = title,
        x = "UMAP 1",
        y = "UMAP 2",
        color = col) +
      theme(legend.position = "bottom") +
      scale_color_manual(values = color.mapping) +
      scale_size_continuous(range = c(1, point.size)) +
      guides(size = "none")
    
  } else {
    clustplot <- data %>% 
      filter(!label %in% c("Other", "Did not meet threshold") & !is.na(label)) %>% 
        ggplot(aes_string(x = umap_x, y = umap_y, color = "label")) +
          geom_point(data = data %>% filter(label %in% c("Other", "Did not meet threshold", NA) | is.na(label)), 
            aes(x = !!sym(umap_x), y = !!sym(umap_y), color = label, size = 0.25),
            alpha = 0.4) +
          geom_point(aes(size = n),
            alpha = 0.8) +
          labs(title = title,
            x = "UMAP 1",
            y = "UMAP 2",
            color = col) +
          theme(legend.position = "bottom") +
          scale_color_manual(values = color.mapping) +
          scale_size_continuous(range = c(1, point.size)) +
          guides(size = "none")
  }

  return(clustplot)
  
}

# Generic function to plot CDR3a and CDR3b UMAPs with the top n label values used to color the UMAP plots.
# Commonly, these labels are clonotype, CDR cluster, T cell type, donor, and predicted epitope. 
plotCDRClusters <- function(data, count.a, col.a, count.b, col.b, top = 20, point.size = 10, levels.a = NULL, levels.b = NULL, method = FALSE) {
  clustplota <- plotCDRCluster(data, count.label = count.a, col = col.a, top, point.size, "CDR3a_UMAP_1", "CDR3a_UMAP_2", "CDR3a Clustering", levels.a, method)
  clustplotb <- plotCDRCluster(data, count.label = count.b, col = col.b, top, point.size, "CDR3b_UMAP_1", "CDR3b_UMAP_2", "CDR3b Clustering", levels.a, method)

  return(list(clustplota, clustplotb))
}

# Wrapper for clonotype overlays
plotClonotypeCluster <- function(data, col = "clonotype", top = 20, point.size = 10) {
  count <- countLabelValues(data, col)
  data <- data %>%
    rename_with(~"label", .cols = all_of(col)) %>%
    left_join(count, by = "label") %>%
    mutate(label = ifelse(label %in% (count %>% top_n(top, n) %>% pull(label)), 
      label, "Other"))
  labels <- unique(data$label)
  levels <- labels[order(as.numeric(gsub("clonotype", "", labels)))] 
  data$label <- factor(data$label, levels = unique(data$label))
  colors <- viridis::turbo(length(labels))
  names(colors) <- labels
  colors <- c(colors, "Other" = "#717070", "Did not meet threshold" = "#717070")
  data$color <- colors[data$label]

  return(plotCDRClusters(data, count.a = count, col.a = col, count.b = count, col.b = col, top, point.size, levels.a = levels, levels.b = levels))

}

# Wrapper for CDR3 cluster overlays
plotCDR3Cluster <- function(data, col.a = "CDR3a_cluster_0.5", col.b = "CDR3b_cluster_0.5", top = 20, point.size = 2) {
  count.a <- countLabelValues(data, "CDR3a_cluster_0.5")
  count.b <- countLabelValues(data, "CDR3b_cluster_0.5")
  data <- data %>%
    mutate(label.a = !!sym(col.a), label.b = !!sym(col.b)) %>%
    select(-!!sym(col.a), -!!sym(col.b)) %>%
    #rename(label.a = !!sym(col.a), label.b = !!sym(col.b)) #%>%
    left_join(count.a, join_by("label.a" == "label"), suffix = c("", ".a")) %>%
    left_join(count.b, join_by("label.b" == "label"), suffix = c(".a", ".b")) %>%
    mutate(label.a = ifelse(label.a %in% (count.a %>% top_n(top, n) %>% pull(label)), 
      label.a, "Other"),
      label.b = ifelse(label.b %in% (count.b %>% top_n(top, n) %>% pull(label)), 
        as.character(label.b), "Other"))
  labels.a <- unique(data$label.a)
  levels.a <- labels.a[order(as.numeric(labels.a))]
  colors.a <- viridis::turbo(length(labels.a))
  names(colors.a) <- labels.a
  colors.a[["Other"]] <- "#717070"
  colors.a <- c(colors.a, "Other" = "#717070", "Does not meet threshold" = "#717070")
  data$color.a <- colors.a[data$label.a]
  labels.b <- unique(data$label.b)
  levels.b <- labels.b[order(as.numeric(labels.b))]
  colors.b <- viridis::turbo(length(labels.b))
  names(colors.b) <- labels.b
  colors.b <- c(colors.b, "Other" = "#717070", "Does not meet threshold" = "#717070")
  data$color.b <- colors.b[data$label.b]

  clustplota <- plotCDRCluster((data %>% 
    mutate(label = label.a, n = n.a) %>% select(-label.a, -n.a) %>% dplyr::rename(color = color.a)), 
    count.label = count.a, col = col.a, top, point.size, "CDR3a_UMAP_1", "CDR3a_UMAP_2", "CDR3a Clustering", levels.a)
  clustplotb <- plotCDRCluster((data %>% 
    mutate(label = label.b, n = n.b) %>% select(-label.b, -n.b) %>% dplyr::rename(color = color.b)),
    count.label = count.b, col = col.b, top, point.size, "CDR3b_UMAP_1", "CDR3b_UMAP_2", "CDR3b Clustering", levels.b)

  return(list(clustplota, clustplotb))
}

# Wrapper for epitope prediction overlays
plotEpitopeCluster <- function(data, col = "sTCR.epitope", data.subset, top = 20, point.size = 10, colors, method = FALSE) {
  count <- countLabelValues(data, col)
  data <- data %>%
    rename_with(~"label", .cols = all_of(col)) %>%
    left_join(count, by = "label") %>%
    mutate(label = ifelse(clone.id %in%
      (data.subset %>% pull(clone.id)), label, "Did not meet threshold")) %>%
    arrange(label == "Did not meet threshold")
  labels <- unique(data$label)
  levels <- c(setdiff(labels, c("Did not meet threshold", "Other")), "Other", "Did not meet threshold")
  data$color <- colors[data$label]
  return(plotCDRClusters(data, count.a = count, col.a = col, count.b = count, col.b = col, top, point.size, levels.a = NULL, levels.b = NULL, method))
}


## Method set for plotting CDR3 motifs ##

# Makes individual motif plots per chain
plotMotif <- function(data, count.label, data.col, cluster.col, chain, levels = NULL) {
  data <- data %>%
    mutate(label = ifelse(label != "Other", paste0(label, "\nn = ", n), "Other"),
      chain = chain) %>%
    filter(label != "Other")

  if (!is.null(levels)) {
    data <- relevelPlot(data, count.label, levels)
  }

  motifplot <- data %>% plot_motifs(data_col = data.col, cluster_col = "label", chain_col = "chain", chain = chain)
  return(motifplot)
  
}

# Generic function that relies on a previous count implementation
plotCDRMotifs <- function(data, count.a, data.col.a = "CDR3a", cluster.col.a, count.b, data.col.b = "CDR3b", cluster.col.b, top = 20, levels.a = NULL, levels.b = NULL) {
  motifplota <- plotCDRMotif(data, count.label = count.a, data.col = data.col.a, cluster.col = cluster.col.a, chain = "A", levels.a)
  motifplotb <- plotCDRMotif(data, count.label = count.b, data.col = data.col.b, cluster.col = cluster.col.b, chain = "B", levels.b)

  return(list(motifplota, motifplotb))
}

# Wrapper for plotting CDR3 clusters with motifs
plotCDRClusterMotifs <- function(data, data.col.a = "CDR3a", cluster.col.a = "CDR3a_cluster_0.5", data.col.b = "CDR3b", cluster.col.b = "CDR3b_cluster_0.5", top = 20) {
  count.a <- countLabelValues(data, cluster.col.a)
  count.b <- countLabelValues(data, cluster.col.b)
  data <- data %>%
    mutate(label.a = !!sym(cluster.col.a), label.b = !!sym(cluster.col.b)) %>%
    select(-!!sym(cluster.col.a), -!!sym(cluster.col.b)) %>%
    #rename("label.a" = cluster.col.a, "label.b" = cluster.col.b) %>%
    left_join(count.a, join_by("label.a" == "label"), suffix = c("", ".a")) %>%
    left_join(count.b, join_by("label.b" == "label"), suffix = c(".a", ".b")) %>%
    mutate(label.a = ifelse(label.a %in% (count.a %>% top_n(top, n) %>% pull(label)), 
      label.a, "Other"),
      label.b = ifelse(label.b %in% (count.b %>% top_n(top, n) %>% pull(label)), 
        label.b, "Other"))
  labels.a <- unique(data$label.a)
  levels.a <- labels.a[order(as.numeric(labels.a))]
  labels.b <- unique(data$label.b)
  levels.b <- labels.b[order(as.numeric(labels.b))]

  motifplota <- plotMotif((data %>% 
    mutate(label = label.a, n = n.a) %>% select(-label.a, -n.a)),
    count.label = count.a, data.col = data.col.a, cluster.col = cluster.col.a, chain = "A", levels.a)
  motifplotb <- plotMotif((data %>% 
    mutate(label = label.b, n = n.b) %>% select(-label.b, -n.b)), 
    count.label = count.b, data.col = data.col.b, cluster.col = cluster.col.b, chain = "A", levels.b)

  return(list(motifplota, motifplotb))
}

# Wrapper for plotting CDR3 motifs per predicted epitope
plotEpitopeMotifs <- function(data, data.col.a = "CDR3a", cluster.col.a = "Epitope", data.col.b = "CDR3b", cluster.col.b = "Epitope", top = 20) {
  count.a <- countLabelValues(data, cluster.col.a)
  count.b <- countLabelValues(data, cluster.col.b)
  data <- data %>%
    mutate(label.a = !!sym(cluster.col.a), label.b = !!sym(cluster.col.b)) %>%
    select(-!!sym(cluster.col.a), -!!sym(cluster.col.b)) %>%
    #rename("label.a" = cluster.col.a, "label.b" = cluster.col.b) %>%
    left_join(count.a, join_by("label.a" == "label"), suffix = c("", ".a")) %>%
    left_join(count.b, join_by("label.b" == "label"), suffix = c(".a", ".b")) %>% 
    filter(n.a > 2, n.b > 2)
  motifplota <- plotMotif((data %>% 
    mutate(label = label.a, n = n.a) %>% select(-label.a, -n.a)),
    count.label = count.a, data.col = data.col.a, cluster.col = cluster.col.a, chain = "A")
  motifplotb <- plotMotif((data %>% 
    mutate(label = label.b, n = n.b) %>% select(-label.b, -n.b)), 
    count.label = count.b, data.col = data.col.b, cluster.col = cluster.col.b, chain = "A")

  return(list(motifplota, motifplotb))

}


## Method set for plotting V/J gene usage ##

# Makes individual V/J Sankey plots
# Data must be filtered to subset prior to calling function
# If data is filtered, this method can be called in standalone version
plotVJSankey <- function(data) {
  vj.freq <- data %>% 
    count(AV, AJ, BV, BJ)
  nodes <- data.frame(name = unique(c(vj.freq$AV, vj.freq$AJ, vj.freq$BV, vj.freq$BJ)))
  nodes$name <- sortAlleles(nodes$name)
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

# General function to plot V/J sankey given a dataframe and column to facilitate 'top n' subsetting
# Common use is to subset by top clonotypes or CDR3 clusters
plotClusterVJSankey <- function(data, col, top = 20) {
  count <- countLabelValues(data, col)
  data <- data %>%
    mutate(label = !!sym(col)) %>%
    select(-!!sym(col)) %>%
    #rename("label" = col) %>%
    filter(label %in% (count %>% top_n(top, n) %>% pull(label)))
  return(plotVJSankey(data))
}


