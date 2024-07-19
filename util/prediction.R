
##########################################
# SwarmTCR training + validation methods #
##########################################


## Basic helper functions ##



## Reference set preparation ##

# There are three files per epitope required to train and validate SwarmTCR: 
# a reference file (distinct from train and test) which contains only positive binding pairs (flag 1),
# and training and test files which contains positive and negative binding pairs (flag 1 and 0, respectively)
splitReferenceData <- function(data, epitopes, directory = "./reference-data/", ref.prop = 0.3, train.prop = 0.7, test.prop = 0.3, n.folds = 25) {

  data <- data %>% 
    formatSwarmTCR(keep_epitopes = TRUE)

  # Create a list of epitope-specific dataframes for each epitope in epitopes
  epitope.dfs <- list()
  for (epitope in epitopes) {
    epitope.dfs[[epitope]] <- data %>% 
      filterEpitope(epitope) %>%
      select(-Epitope)
  }

  # Split the epitope reference datasets into train, test, and ref sets
  for (epitope in epitopes) {
    epitope.df <- epitope.dfs[[epitope]] 

    for (i in 1:n.folds) {
      ref.data <- epitope.df %>% 
        slice_sample(n = floor(ref.prop * nrow(epitope.df))) %>% 
        mutate(Flag = 1) %>%
        select(TCR_ID, Flag, CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b)
      train.test.data <- data %>%
        mutate(Flag = ifelse(str_detect(Epitope, epitope), 1, 0)) %>%
        select(-Epitope) %>%
        anti_join(ref.data, by = "TCR_ID") %>%
        select(TCR_ID, Flag, CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b)
      train.data <- train.test.data %>% 
        slice(createDataPartition(train.test.data$Flag, p = train.prop, list = FALSE))
      test.data <- train.test.data %>%
        anti_join(train.data, by = "TCR_ID")
      
      dir.create(paste0(directory, epitope, "/fold", i, "/"), recursive = TRUE, showWarnings = FALSE)
      write_delim(ref.data, paste0(directory, epitope, "/fold", i, "/", "Reference.txt"))
      write_delim(train.data, paste0(directory, epitope, "/fold", i, "/", "Train.txt"))
      write_delim(test.data, paste0(directory, epitope, "/fold", i, "/", "Test.txt"))
    }
  }
}


# To use SwarmTCRClassify after training, each epitope needs a reference set containing all positive binding pairs (no split)
# Or, use the pre-trained weights and references provided in SwarmClassify
makeClassifierReference <- function(data, epitopes, directory = "./reference-data/") {
  for (epitope in epitopes) {
    epitope.df <- data %>% 
      filterEpitope(epitope) %>%
      select(-Epitope) %>% 
      formatSwarmTCR(flag = 1, keep_epitopes = FALSE)
    dir.create(paste0(directory, epitope, "/"), recursive = TRUE, showWarnings = FALSE)
    write_delim(epitope.df, paste0(directory, epitope, "/", epitope, "_Classifier-Reference.txt"))
  }
}

#######################################
# SwarmTCR results processing methods #
#######################################

## Functions to process each component of SwarmTCR results ##

# Process individual refWeights files
readRefWeights <- function(refweight.file) {
  lines <- readLines(refweight.file)
  filtered.lines <- grep("^#\\d+", lines, value = TRUE)
  refweights <- data.frame()
  final.line <- filtered.lines[length(filtered.lines)]
  # Split the line into components
  components <- strsplit(final.line, " ")[[1]]
  # Remove the first component ('#digit') and the '--'
  components <- components[-c(1, 3)]
  weights <- as.numeric(components)
  weights <- data.frame(
  Score = weights[1],
  CDR1a = weights[2],
  CDR2a = weights[3],
  CDR2_5a = weights[4],
  CDR3a = weights[5],
  CDR1b = weights[6],
  CDR2b = weights[7],
  CDR2_5b = weights[8],
  CDR3b = weights[9])
  refweights <- rbind(refweights, weights)
  
  return(refweights)
}

# Process all refWeights file
# Returns lists of per-epitope reference weights dataframes 
processRefWeights <- function(directory, epitopes) {
  refweights <- list()
  for (epitope in epitopes) {
    epitope.directory <- paste0(directory, epitope, "/")
    refweight.files <- list.files(epitope.directory, pattern = "refWeights", full.names = TRUE, recursive = TRUE)
    epitope.refweights <- data.frame()
    for (refweight.file in refweight.files) {
      refweight <- readRefWeights(refweight.file)
      epitope.refweights <- rbind(epitope.refweights, refweight)
    }
    refweights[[epitope]] <- epitope.refweights
  }
  return(refweights)
}

# Process individual pr files
readPR <- function(pr.file, run) {
  pr <- read.table(pr.file, sep = "\t", check.names = FALSE, header = TRUE) %>%
    mutate(Run = run)
  return(pr)
}

# Process all pr files
# Returns two lists of precision/recall dataframes (swarmTCR + TCRdist)
processPR <- function(directory, epitopes) {
  pr.swarmTCR <- list()
  pr.TCRdist <- list()
  for (epitope in epitopes) {
    epitope.directory <- paste0(directory, epitope, "/")
    pr.files.swarmTCR <- list.files(epitope.directory, pattern  = "val_SwarmTCR_out\\.txt$", full.names = TRUE, recursive = TRUE)
    pr.files.TCRdist <- list.files(epitope.directory, pattern  = "val_TCRdist_out\\.txt$", full.names = TRUE, recursive = TRUE)
    epitope.pr.swarmTCR <- data.frame()
    epitope.pr.TCRdist <- data.frame()
    for (i in 1:length(pr.files.swarmTCR)) {
      pr <- readPR(pr.files.swarmTCR[i], i)
      epitope.pr.swarmTCR <- rbind(epitope.pr.swarmTCR, pr)
      pr <- readPR(pr.files.TCRdist[i], i)
      epitope.pr.TCRdist <- rbind(epitope.pr.TCRdist, pr)
    }
    pr.swarmTCR[[epitope]] <- epitope.pr.swarmTCR
    pr.TCRdist[[epitope]] <- epitope.pr.TCRdist
  }
  return(list(pr.swarmTCR, pr.TCRdist))
}

# Process individual scores files
readScores <- function(scores.file, run) {
  scores <- read.table(scores.file, sep = "\t", check.names = FALSE, header = FALSE) %>%
    mutate(Run = run)
  return(scores)
}
# Process all scores files
# Returns list of validation scores dataframes
processValidationScores <- function(directory, epitopes) {
  scores <- list()
  for (epitope in epitopes) {
    epitope.directory <- paste0(directory, epitope, "/")
    scores.files <- list.files(epitope.directory, pattern = "\\.txt\\.scores$", full.names = TRUE, recursive = TRUE)
    epitope.scores <- data.frame()
    for (i in 1:length(scores.files)) {
      score <- readScores(scores.files[i], i)
      epitope.scores <- rbind(epitope.scores, score)
    }
    colnames(epitope.scores) <- c("TCR_ID", "Flag", "Score", "Run")
    scores[[epitope]] <- epitope.scores
  }
  return(scores)
}

## Wrapper function to process all results at once ##
processSwarmTCRResults <- function(directory, epitopes) {
  results <- list()
  results[["refweights"]] <- processRefWeights(directory, epitopes)
  pr <- processPR(directory, epitopes)
  results[["pr.swarmTCR"]] <- pr[[1]]
  results[["pr.TCRdist"]] <- pr[[2]]
  results[["scores"]] <- processValidationScores(directory, epitopes)
  return(results)
}

## Functions to set SwarmClassify parameters ##

# Function set for setting reference weights
setRefWeights <- function(refweights.dfs, epitopes.df, directory = "./reference-data/", append.directory = "../reference-data/", method = "mediod", threshold = NULL) {
  set.weights <- data.frame()
  all.weights <- data.frame()

  epitopes <- epitopes.df %>% pull(Epitope)

  # Get per-epitope best and mediod reference weights
  for (epitope in epitopes) {
    refweights <- refweights.dfs[[epitope]]
    best.weights <- getBestWeights(refweights, epitope) %>% 
        slice(1)
    mediod.weights <- getMediodWeights(refweights, epitope) %>% 
        slice(1)
    if (method == "best") {
      set.weights <- rbind(set.weights, best.weights)
    } else if (method == "mediod") {
      set.weights <- rbind(set.weights, mediod.weights)
    }
    all.weights <- rbind(all.weights, (best.weights %>% 
      mutate(Method = "best")) %>%
      rbind(mediod.weights %>%
        mutate(Method = "mediod")))
  }

  # Filter reference weights by threshold for refWeights file and final epitopes file, if given
  if (!is.null(threshold)) {
    set.weights <- set.weights %>% 
      filter(Score >= threshold)
    epitopes.df <- epitopes.df %>%
      filter(Epitope %in% (set.weights %>% pull(Epitope)))
    epitopes <- epitopes.df %>% pull(Epitope) 

    # Write final epitopes to file if threshold is given
    epitope.file <- list.files(directory, pattern = "epitopes\\.csv$", full.names = FALSE, recursive = FALSE)
    write.csv(epitopes.df, paste0(directory, "final-", epitope.file), row.names = FALSE)
  }

  # Write final refWeights file
  makeRefWeightsFile(set.weights, directory, append.directory)

  # Return reference weights results dataframes list (set refWeights, all refWeights, valid epitopes)
  results <- list()
  results[["set.weights"]] <- set.weights
  results[["all.weights"]] <- all.weights
  results[["final.epitopes"]] <- epitopes.df
  return(results)
}

getBestWeights <- function(refweights, epitope) {
  setweight <- refweights %>% 
    filter(Score == max(Score)) %>%
    mutate(Epitope = epitope)

  return(setweight)
}

getMediodWeights <- function(refweights, epitope) {
  cdr.columns <- grep("CDR", names(refweights), value = TRUE)
  weightdist <- refweights %>% 
    rowwise() %>%
    mutate(across(all_of(cdr.columns), ~ {
      column.values <- refweights[[cur_column()]]
      dists <- sqrt((column.values - .)^2)
      sum(dists)
    }, .names = "{.col}.dist")) %>%
    ungroup() %>%
    mutate(total.dist = rowSums(select(., ends_with(".dist"))),
           Epitope = epitope) 

  setweight <- weightdist %>%
    filter(total.dist == min(total.dist)) %>%
    select(-ends_with(".dist"))
  
  return(setweight)
}

makeRefWeightsFile <- function(set.weights, directory, append.directory = "../reference-data/") {
  set.weights <- set.weights %>% 
    mutate(FILE_NAME = paste0(append.directory, Epitope, "/", Epitope, "_Classifier-Reference.txt")) %>%
    select(FILE_NAME, Epitope, CDR1a, CDR2a, CDR2_5a, CDR3a, CDR1b, CDR2b, CDR2_5b, CDR3b)
  write_delim(set.weights, paste0(directory, "refWeights.txt"))
}

# Function set for setting score thresholds
setScoreThresholds <- function(pr.dfs, scores.dfs, epitopes, directory = "./reference-data/", method = "precision") {
  thresholds <- data.frame()
  for (epitope in epitopes) {
    pr <- pr.dfs[[epitope]]
    scores <- scores.dfs[[epitope]]
    if (method == "precision") {
      threshold <- getScorePrecision(pr, scores, epitope)
    } else if (method == "f1") {
      threshold <- getScoreF1(pr, scores, epitope)
    }
    thresholds <- rbind(thresholds, threshold)
  }
  write.csv(thresholds, paste0(directory, "scoreThresholds.csv"), row.names = FALSE)
  return(thresholds)
}


getScoreF1 <- function(pr, scores, epitope) {
  optimal <- pr %>%
    mutate(F1 = 2 * Precision * Recall / (Precision + Recall)) %>% 
    filter(F1 == max(F1)) 
  f1 <- optimal %>%
    pull(F1)
  precision <- optimal %>%
    pull(Precision)
  recall <- optimal %>%
    pull(Recall)
  scores <- scores %>%
    filter(Flag == 1) %>% 
    arrange(desc(Score)) 
  threshold <- scores %>%
    slice(round(recall * nrow(scores))) %>%
    pull(Score)
  return(data.frame(Epitope = epitope, F1 = f1, Precision = precision, Recall = recall, Threshold = threshold))
}

getScorePrecision <- function(pr, scores, epitope) {
  optimal <- pr %>%
    mutate(F1 = 2 * Precision * Recall / (Precision + Recall)) %>%
    arrange(desc(Precision), desc(Recall)) %>%
    slice_head(n = 1)
  f1 <- optimal %>%
    pull(F1)
  precision <- optimal %>%
    pull(Precision)
  recall <- optimal %>%
    pull(Recall)
  scores <- scores %>%
    filter(Flag == 1) %>% 
    arrange(desc(Score)) 
  threshold <- scores %>%
    slice(round(recall * nrow(scores))) %>%
    pull(Score)
  return(data.frame(Epitope = epitope, F1 = f1, Precision = precision, Recall = recall, Threshold = threshold))
}




############################################
# SwarmClassify results processing methods #
############################################

processSwarmClassify <- function(results.path, exp.name, cdrs.path = NULL) {
      exp.results <- import(paste0(results.path, exp.name, "_CDRs.txt.results"))
    if (!is.null(input.path)) {
        exp.cdrs <- import(paste0(input.path, exp.name, "_CDRs.txt"))
        exp.results <- merge(exp.cdrs, exp.results, by = "TCR_ID")
    }
    return(exp.results)
}
