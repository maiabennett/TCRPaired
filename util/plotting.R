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

## General functions to assist in plotting
# Counting the number of times a labeling (overlay) value appears in the data to facilitate top n + "other" labeling
countLabelValues <- function(data, col) {
  count.label <- data %>%
    count(!!sym(col)) %>%
    rename("label" = col)
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
plotCDRCluster <- function(data, count.label, col, top, point.size, umap_x, umap_y, title, levels = NULL) {
  data <- data %>%
    mutate(label = ifelse(label != "Other", paste0(label, "\nn = ", n), "Other"))

  if (!is.null(levels)) {
    data <- relevelPlot(data, count.label, levels)
  }

  data %>%
    ggplot(aes_string(x = umap_x, y = umap_y, color = "label")) +
    geom_point(aes(size = n), alpha = 0.7) +
    labs(title = title,
      x = "UMAP 1",
      y = "UMAP 2",
      color = col) +
    theme(legend.position = "bottom") +
    scale_color_viridis(option = "F", discrete = TRUE, begin = 1, end = 0.2) +
    scale_size_continuous(range = c(1, point.size))
}

# Generic function to plot CDR3a and CDR3b UMAPs with the top n label values used to color the UMAP plots.
# Commonly, these labels are clonotype, CDR cluster, T cell type, donor, and predicted epitope. 
plotCDRClusters <- function(data, count.a, col.a, count.b, col.b, top = 20, point.size = 10, levels.a = NULL, levels.b = NULL) {
  clustplota <- plotCDRCluster(data, count.label = count.a, col = col.a, top, point.size, "CDR3a_UMAP_1", "CDR3a_UMAP_2", "CDR3a Clustering", levels.a)
  clustplotb <- plotCDRCluster(data, count.label = count.b, col = col.b, top, point.size, "CDR3b_UMAP_1", "CDR3b_UMAP_2", "CDR3b Clustering", levels.a)

  return(list(clustplota, clustplotb))
}

# Wrapper for clonotype overlays
plotClonotypeCluster <- function(data, col = "clonotype", top = 20, point.size = 10) {
  count <- countLabelValues(data, col)
  data <- data %>%
    rename("label" = col) %>%
    left_join(count, by = "label") %>%
    mutate(label = ifelse(label %in% (count %>% top_n(top, n) %>% pull(label)), 
      label, "Other"))
  labels <- unique(data$label)
  levels <- labels[order(as.numeric(gsub("clonotype", "", labels)))] 

  return(plotCDRClusters(data, count.a = count, col.a = col, count.b = count, col.b = col, top, point.size, levels.a = levels, levels.b = levels))

}

# Wrapper for CDR3 cluster overlays
plotCDR3Cluster <- function(data, col.a = "CDR3a_cluster_0.5", col.b = "CDR3b_cluster_0.5", top = 20, point.size = 2) {
  count.a <- countLabelValues(data, col.a)
  count.b <- countLabelValues(data, col.b)
  data <- data %>%
    rename("label.a" = col.a, "label.b" = col.b) %>%
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

  clustplota <- plotCDRCluster((data %>% rename("label" = "label.a", "n" = "n.a")), count.label = count.a, col = col.a, top, point.size, "CDR3a_UMAP_1", "CDR3a_UMAP_2", "CDR3a Clustering", levels.a)
  clustplotb <- plotCDRCluster(data %>% rename("label" = "label.b", "n" = "n.b"), count.label = count.b, col = col.b, top, point.size, "CDR3b_UMAP_1", "CDR3b_UMAP_2", "CDR3b Clustering", levels.b)

  return(list(clustplota, clustplotb))
}

# Wrapper for epitope prediction overlays
#plotEpitopeCluster <- function()


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
    rename("label.a" = cluster.col.a, "label.b" = cluster.col.b) %>%
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

  motifplota <- plotMotif((data %>% rename("label" = "label.a", "n" = "n.a")), count.label = count.a, data.col = data.col.a, cluster.col = cluster.col.a, chain = "A", levels.a)
  motifplotb <- plotMotif((data %>% rename("label" = "label.b", "n" = "n.b")), count.label = count.b, data.col = data.col.b, cluster.col = cluster.col.b, chain = "A", levels.b)

  return(list(motifplota, motifplotb))
}

# Wrapper for plotting CDR3 motifs per predicted epitope
#plotEpitopeMotifs <- function()


## Method set for plotting V/J gene usage ##

# Makes individual V/J Sankey plots
# Data must be filtered to subset prior to calling function
# If data is filtered, this method can be called in standalone version
plotVJSankey <- function(data) {
  vj.freq <- data %>% 
    count(AV, AJ, BV, BJ)
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
    rename("gene" = "name")

  sankeyplot <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source", 
    Target = "target", Value = "value", NodeID = "name", sinksRight = FALSE, LinkGroup = "gene")
}

# General function to plot V/J sankey given a dataframe and column to facilitate 'top n' subsetting
# Common use is to subset by top clonotypes or CDR3 clusters
plotClusterVJSankey <- function(data, col, top = 20) {
  count <- countLabelValues(data, col)
  data <- data %>%
    rename("label" = col) %>%
    filter(label %in% (count %>% top_n(top, n) %>% pull(label)))
  return(plotVJSankey(data))
}

# General function to plot V/J sankey given a dataframe, a column, and a precise target value to subset by 
# Common use is to subset by epitope, epitope species