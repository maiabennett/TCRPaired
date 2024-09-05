


#####################
# Filtering methods #
#####################

## Filtering entries by distinct sets of sequences to determine unique receptors ##

# Filter to unique full sequences
filterFullSeq <- function(data) {
  data <- data %>% 
    distinct(Epitope, full.seq, .keep_all = TRUE) %>%
    drop_na(full.seq) %>%
    filter(full.seq != "")
}

# Filter to unique CDR sequences
filterCDRSeq <- function(data) {
    if (any(names(data) == "Epitope")) {
      data <- data %>% 
        distinct(CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b, Epitope, .keep_all = TRUE) %>%
        drop_na(CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b) %>%
        filter_at(vars(CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b), 
          all_vars(. != ""))
    } else {
      data <- data %>% 
        distinct(CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b, .keep_all = TRUE) %>%
        drop_na(CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b) %>%
        filter_at(vars(CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b), 
          all_vars(. != ""))
    }
}

# Filter to unique V/J + CDR3 entries
filterVJCDR3 <- function(data) {
  if (any(names(data) == "Epitope")) {
    data <- data %>% 
      distinct(AV, CDR3a, AJ, BV, CDR3b, BJ, Epitope, .keep_all = TRUE) %>% 
      drop_na(AV, CDR3a, AJ, BV, CDR3b, BJ, Epitope) %>% 
      filter_at(vars(AV, CDR3a, AJ, BV, CDR3b, BJ, Epitope), 
        all_vars(. != ""))
  } else {
    data <- data %>% 
      distinct(AV, CDR3a, AJ, BV, CDR3b, BJ, .keep_all = TRUE) %>% 
      drop_na(AV, CDR3a, AJ, BV, CDR3b, BJ) %>% 
      filter_at(vars(AV, CDR3a, AJ, BV, CDR3b, BJ), 
        all_vars(. != ""))
  }
}


## Filtering by common column values ##

# Filter to high-confidence entries 
# Default is VDJdb score > 0 + any clone.id identifiers passed (e.g., crossreactive)
filterConfidence <- function(data, score = 1, ...) {
  sources <- c(...)
  data <- data %>%
    filter(Score >= score | str_detect(clone.id, paste0("^(", paste(sources, collapse = "|"), ")")))
}

# Filter by epitope gene
filterEpitopeGene <- function(data, epitope.gene, epitope.species = NULL) {
  if (is.null(epitope.species)) {
    data <- data %>% 
      filter(str_detect(Epitope.gene, epitope.gene))
  } else {
    data <- data %>% 
      filter(str_detect(Epitope.gene, epitope.gene) & str_detect(Epitope.species, epitope.species))
  }
}

# Filter by epitope species
filterEpitopeSpecies <- function(data, epitope.species) {
  data <- data %>%
    filter(str_detect(Epitope.species, epitope.species))
}


## Filtering by similarity and distance ##
## Internal methods ##

# computeSimilarity function calculates the similarity between each pair of items in the data based on their 'column' values.
# It uses the Levenshtein distance (edit distance) to compute a similarity matrix, then transforms this into similarity levels.
# The function adds two new columns to the data: 'similarity.level' and 'max.similarity.level', where 'max.similarity.level'
# represents the maximum similarity of each item with any other item.
# Compare.x will always be a column of data to compare
# Compare.y can be a column or a sequence to compare against
computeSimilarity <- function(data, compare.x, compare.y, similarity = 0.9) {
  #column <- data %>% select(-clone.id) %>% pull()
    data <- data %>% 
        mutate(dist.matrix = stringdistmatrix(compare.x, compare.y, method = "lv"))
    if(is.data.frame(compare.y)) {
      diag(data$dist.matrix) <- NA
    }
    data <- data %>%
        mutate(similarity.level = 1 - dist.matrix / nchar(compare.x),
            max.similarity.level = apply(similarity.level, 1, max, na.rm = TRUE))
    return(data)
}

# computeRepresentatives function identifies clusters of items in the data that are similar to each other above a specified threshold.
# It first filters items based on their 'max.similarity.level', then computes a similarity matrix for the filtered items.
# This matrix is used to perform hierarchical clustering, and clusters are defined by cutting the dendrogram at a height that reflects
# the similarity threshold. This process groups items into clusters based on their similarity, identifying representative items or groups.
computeRepresentatives <- function(data, similarity = 0.9) {
    clusters <- unique(data$cluster_idx)

    rep.data <- data.frame()

    for (cluster in clusters) {
        cluster.data <- data %>%
            filter(cluster_idx == cluster) 
        cluster.column <- cluster.data %>% 
            select(-cluster_idx, -n_cluster) %>%
            pull()
        dist.matrix <- stringdistmatrix(cluster.column, cluster.column, method = "lv")
        # Use mediod as representative
        sum.distances <- apply(dist.matrix, 1, sum)
        medoid <- cluster.data[which.min(sum.distances),]
        rep.data <- rbind(rep.data, medoid)
    }

    return(rep.data)

  #sim <- data %>% 
  #  filter(max.similarity.level >= similarity) 
  #if (nrow(sim) < 2) {
  #  return(sim)
  #} else {
  #  column <- sim %>%
  #    select(-dist.matrix, -similarity.level, -max.similarity.level, -clone.id) %>% pull()
  #  dist.matrix <- stringdistmatrix(column, column, method = "lv")
  #  sim.matrix <- 1 - dist.matrix / max(dist.matrix)
  #  sim.clust <- hclust(as.dist(1-sim.matrix))
  #  sim.clusters <- cutree(sim.clust, h = 1-similarity)
  #  sim <- sim %>%
  #    mutate(cluster = sim.clusters) %>%
  #    group_by(cluster) %>%
  #    arrange(desc(max.similarity.level)) %>%
  #    slice(1) %>%
  #    ungroup()

  #  return(sim %>% select(-cluster))
  #}
}

computeDistance <- function(data, compare.x, compare.y, method = "lv") {
    data <- data %>% 
        mutate(dist = stringdistmatrix(compare.x, compare.y, method = method))
    return(data)
}

## Wrapper methods ##

# Filter by degree of similarity: called on its own or by other methods (epitope, motif)
# When keep is true, filter out highly similar sequences while keeping one representative for each set of sequences which are 90% similar to one another
# When keep and above are false, filter out all rows which have highly similar sequences
# When above is true, filter out all rows which have sequences that are less than 90% similar to one another, 
# essentially performing a motif search but with percent similarity rather than exact allowable distance

# Re-working based on CD-HIT and CellaRepetorium
# With the idea that we can lose some sequences with single-AA missing values (X, #, etc.) for sake of finding similar sequences in large datasets
# As such, there's two wrappers from the previous one: 
# filterSequenceSimilarity and filterMotifSimilarity
# The former is for filtering by full sequences, producing massively polynomial matrices, the latter is for filtering by motifs (CDR3a, CDR3b, Epitopes, etc.), producing fixed-size matrices

filterSequenceSimilarity <- function(data, column, similarity = 0.9) {
    filter.by <- data %>%
        checkAASequence() %>% 
        select(clone.id, column)
    filter.values <- filter.by %>% 
        select(-clone.id) %>%
        pull() 
    names(filter.values) <- filter.by$clone.id
    sim.clusters <- filter.values %>%
        AAStringSet() %>%
        cdhit(identity = similarity)
    
    # Store dissimilar seqs below threshold
    dissim.seqs <- sim.clusters %>%
        filter(n_cluster == 1)
    
    # Compute representatives for clusters with > 1 sequence at similarity threshold
    sim.seqs <- sim.clusters %>%
        filter(n_cluster > 1) 
    rep.seqs <- sim.seqs %>% 
        computeRepresentatives(similarity)

    # Combine dissimilar and representative sequences
    seqs <- rbind(dissim.seqs, rep.seqs)

    # Join back to original data where clone.id = query_name
    data <- data %>%
        inner_join(seqs, by = c("clone.id" = "query_name")) %>%
        select(-cluster_idx, -n_cluster, -seq)
    return(data)    
}

filterMotifSimilarity <- function(data, column, motif, similarity = 0.9, above = TRUE) {
    filter.by <- data %>% select(clone.id, column)
    #if (is.null(motif)) {
    #    sim <- computeSimilarity(filter.by, 
    #        (filter.by %>% select(-clone.id) %>% pull()), 
    #        (filter.by %>% select(-clone.id) %>% pull()), 
    #        similarity)
    #} else {
    sim <- computeSimilarity(filter.by, (filter.by %>% select(-clone.id) %>% na.omit() %>% pull()), motif, similarity)
    #}
    #if (representative) {
    #    reps <- computeRepresentatives(sim, similarity)
    #    sim <- rbind(reps, sim %>% filter(max.similarity.level < similarity))
    #}
    if (above) {
      sim <- sim %>% 
        filter(max.similarity.level >= similarity)
    } else {
      sim <- sim %>% 
        filter(max.similarity.level < similarity)
    }
    data <- data %>% inner_join(sim %>% select(-dist.matrix, -similarity.level, -max.similarity.level, clone.id), by = "clone.id")

    return(data) 
}

filterDistance <- function(data, column, motif = NULL, distance = 0, above = FALSE) {
    filter.by <- data %>% select(clone.id, column)
    if (is.null(motif)) {
        sim <- computeDistance(filter.by, 
        (filter.by %>% select(-clone.id) %>% pull()), 
        (filter.by %>% select(-clone.id) %>% pull()))
    } else {
        sim <- computeDistance(filter.by, (filter.by %>% select(-clone.id) %>% pull()), motif)
    } 
    if (above) {
        sim <- sim %>% 
        filter(dist >= distance)
    } else {
        sim <- sim %>% 
        filter(dist <= distance)
    }
    data <- data %>% inner_join(sim %>% select(-dist, clone.id), by = "clone.id")

    return(data)
}

# Filter to epitope-specific entries
# Default is to filter by exact epitope match
# If distance is provided, filter by allowable sequence distance
# If similarity is provided, filter by allowable sequence similarity
filterEpitope <- function(data, epitope, distance = NULL, similarity = NULL, representative = FALSE, above = TRUE) {
  if (is.null(distance) && is.null(similarity)) {
    data <- data %>%
      filter(str_detect(Epitope, epitope))
  } else if (!is.null(distance)) {
    data <- data %>%
      filterDistance(Epitope, motif = epitope, distance, above = above)
  } else if (!is.null(similarity)) {
    data <- data %>%
      filterMotifSimilarity(Epitope, motif = epitope, similarity, representative = representative, above = above)
  }
}

# Filter to CDR3 motif-specific entries
# Default is to filter by exact motif match
# If distance is provided, filter by allowable sequence distance
# If similarity is provided, filter by allowable sequence similarity
filterCDR3a <- function(data, motif, distance = NULL, similarity = NULL, representative = FALSE, above = TRUE) {
  if (is.null(distance) && is.null(similarity)) {
    data <- data %>%
      filter(str_detect(CDR3a, motif))
  } else if (!is.null(distance)) {
    data <- data %>%
      filterDistance(CDR3a, motif = motif, distance, above = above)
  } else if (!is.null(similarity)) {
    data <- data %>%
      filterMotifSimilarity(CDR3a, motif = motif, similarity, representative = representative, above = above)
  }
  return(data)
}
filterCDR3b <- function(data, motif, distance = NULL, similarity = NULL, representative = FALSE, above = TRUE) {
  if (is.null(distance) && is.null(similarity)) {
    data <- data %>%
      filter(str_detect(CDR3b, motif))
  } else if (!is.null(distance)) {
    data <- data %>%
      filterDistance(CDR3b, motif = motif, distance, above = above)
  } else if (!is.null(similarity)) {
    data <- data %>%
      filterMotifSimilarity(CDR3b, motif = motif, similarity, representative = representative, above = above)
  }
  return(data)
}



