
##' @title Classification forest class
##' @description Subclass for classification forest. 
##' Contains all fields and methods used special for classification forests.
ForestClassification <- setRefClass("ForestClassification", 
  contains = "Forest", 
  fields = list(
    response_levels = "character"), 
  methods = list(
    
    grow = function(num_threads) {      
      treetype <<- "Classification"
      
      ## Create trees
      if (splitobject == "module") {
        trees <<- replicate(num_trees, TreeVarClustersClassification$new())
      } else {
        trees <<- replicate(num_trees, TreeClassification$new())
      }
      
      ## Call parent method
      callSuper(num_threads)
    },
    
    aggregatePredictions = function(predictions) {
      ## For all samples take majority vote over all trees
      result <- apply(predictions, 1, function(x) {
        class_counts <- table(x)
        return(as.numeric(names(class_counts)[which.max.random(class_counts)]))
      })
      
      ## Return as factor
      return(factor(result, levels = 1:length(response_levels), labels = response_levels))
    }, 
    
    predictionErrorForest = function() {
      ## For each tree loop over OOB samples and count classes
      tree_predictions <- sapply(trees, function(x) {
        oob_samples <- x$oob_sampleIDs
        result <- rep(NA, data$nrow)
        result[oob_samples] <- x$predictOOB()
        return(result)
      })

      # print(tree_predictions)
      
      ## Compute majority vote for each sample
      sample_predictions <- apply(tree_predictions, 1, function(x) {
        if (sum(!is.na(x)) > 0) {
          class_counts <- table(x)
          return(as.numeric(names(class_counts)[which.max.random(class_counts)]))
        } else {
          return(NA)
        }
      })

      # print(sample_predictions)
      # print(as.numeric(data$column(1)))
      
      ## Compare predictions with true data
      return(sum(sample_predictions != as.numeric(data$column(1)), na.rm = TRUE) / data$nrow)
    },
    
    predictionErrorTrees = function() {
      ## For each tree loop over OOB samples and count classes
      tree_predictions <- lapply(trees, function(x) {
        oob_samples <- x$oob_sampleIDs
        result <- rep(NA, data$nrow)
        result[oob_samples] <- x$predictOOB()
        return(result)
      })

      ## Calculate OOB prediction errors for all trees
      errors <- sapply(
        1:length(trees),
        function(treeID) {
          return(sum(tree_predictions[[treeID]] != as.numeric(data$column(1)), na.rm = TRUE) / data$nrow)
        }
      )

      ## Return mean error
      return(mean(errors))
    },
    
    predictionErrorForestAndTrees = function() {
      ## For each tree loop over OOB samples and count classes
      tree_predictions <- sapply(trees, function(x) {
        oob_samples <- x$oob_sampleIDs
        result <- rep(NA, data$nrow)
        result[oob_samples] <- x$predictOOB()
        return(result)
      })

      ## Compute majority vote for each sample
      sample_predictions <- apply(tree_predictions, 1, function(x) {
        if (sum(!is.na(x)) > 0) {
          class_counts <- table(x)
          return(as.numeric(names(class_counts)[which.max.random(class_counts)]))
        } else {
          return(NA)
        }
      })

      ## Calculate OOB prediction errors for all trees
      errors_trees <- sapply(
        1:length(trees),
        function(treeID) {
          return(sum(tree_predictions[,treeID] != as.numeric(data$column(1)), na.rm = TRUE) / length(trees[[treeID]]$oob_sampleIDs))
        }
      )

      ## Return errors
      return(list(
        error_forest = sum(sample_predictions != as.numeric(data$column(1)), na.rm = TRUE) / sum(!is.na(sample_predictions)),
        error_trees = mean(errors_trees))
      )
    })
)
