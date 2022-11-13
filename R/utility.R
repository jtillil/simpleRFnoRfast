##' Find index of maximum, breaking ties at random and ignoring NAs.
##' @title Index of maximum with random ties
##' @param x Input vector.
##' @return Index of the maximum value.
##' @author Marvin N. Wright
which.max.random <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  which(rank(x, ties.method = "random", na.last = FALSE) == length(x))
}

##' Find index of minimum, breaking ties at random and ignoring NAs.
##' @title Index of minimum with random ties
##' @param x Input vector.
##' @return Index of the minimum value.
##' @author Marvin N. Wright
which.min.random <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  which(rank(x, ties.method = "random", na.last = TRUE) == 1)
}

##' Order factor levels with correlation approach
##' @title Order factor levels with correlation approach
##' @param y Response factor.
##' @param x Covariate factor.
##' @return Ordered factor levels
##' @author Marvin N. Wright
cor.order <- function(y, x) {
  ## Create contingency table of the nominal outcome with the nominal covariate
  tab <- table(droplevels(y), droplevels(x))
  
  ## Compute correlation matrix of the contingency table (correlation of the covariate levels w.r.t outcome)
  cr <- suppressWarnings(cor(tab))
  cr[is.na(cr)] <- 0                      
  diag(cr) <- NA
  
  ## Start with a random level and select as next level the level with highest correlation to the current level (excluding already selected levels)
  num_levels <- nlevels(droplevels(x))
  next_level <- sample(num_levels, 1)
  res <- c(next_level, rep(NA, num_levels - 1))
  for (i in 2:num_levels) {
    cr[, next_level] <- NA
    next_level <- which.max.random(cr[next_level, ])
    res[i] <- next_level
  }
  
  ## Return ordered factor levels
  as.character(levels(droplevels(x))[res])
}

##' Order factor levels with PCA approach
##' @title Order factor levels with PCA approach
##' @param y Response factor.
##' @param x Covariate factor.
##' @return Ordered factor levels
##' @author Marvin N. Wright
##' @references Coppersmith, D., Hong, S.J. & Hosking, J.R. (1999) Partitioning Nominal Attributes in Decision Trees. Data Min Knowl Discov 3:197. \url{https://doi.org/10.1023/A:1009869804967}.
pc.order <- function(y, x) {
  if (nlevels(droplevels(x)) < 2) {
    return(as.character(levels(droplevels(x))))
  }
  
  ## Create contingency table of the nominal outcome with the nominal covariate
  N <- table(droplevels(x), droplevels(y))
  
  ## PCA of weighted covariance matrix of class probabilites
  P <- N/rowSums(N)
  S <- cov.wt(P, wt = rowSums(N))$cov
  pc1 <- prcomp(S, rank. = 1)$rotation
  score <- P %*% pc1
  
  ## Return ordered factor levels
  as.character(levels(droplevels(x))[order(score)])
}

##' Reorder factor columns. Use mean for continuous response, class counts for factors and mean survival for survival response.
##' @title Reorder factor columns
##' @param data Data with factor columns.
##' @return Data with reordered factor columns.
##' @importFrom coin logrank_trafo
##' @author Marvin N. Wright
reorder.factor.columns <- function(data) {
  ## Recode characters and unordered factors
  character.idx <- sapply(data[, -1], is.character)
  ordered.idx <- sapply(data[, -1], is.ordered)
  factor.idx <- sapply(data[, -1], is.factor)
  recode.idx <- character.idx | (factor.idx & !ordered.idx)
  
  ## Numeric response
  response <- data[, 1]
  if (is.factor(response)) {
    num.response <- as.numeric(response)
  } else if ("Surv" %in% class(response)) {
    num.response <- coin::logrank_trafo(response, ties.method = "Hothorn-Lausen")
  } else {
    num.response <- response
  }
  
  ## Recode each column
  data[, -1][, recode.idx] <- lapply(data[, -1][, recode.idx, drop = FALSE], function(x) {
    if (is.factor(response) & nlevels(response) > 2) {
      levels.ordered <- pc.order(y = response, x = x)
    } else {
      ## Order factor levels by num.response
      means <- sapply(levels(x), function(y) {
        mean(num.response[x == y])
      })
      levels.ordered <- as.character(levels(x)[order(means)])
    }
    
    ## Return reordered factor
    factor(x, levels = levels.ordered, ordered = TRUE)
  })
  
  ## Return data
  data
}

##' Convert number to bit vector.
##' @title Number to bit vector
##' @param x Input number.
##' @param length Length of output bit vector.
##' @return Bit vector.
##' @author Marvin N. Wright
as.bitvect <- function(x, length = 32) {
  i <- 1
  string <- numeric(length)
  while(x > 0) {
    string[i] <- x %% 2
    x <- x %/% 2
    i <- i + 1 
  }
  as.logical(string)
}

##' Compute Gini impurity of linear split.
##' @title Gini impurity of linear split
##' @param dat Data frame of predictors.
##' @param label Vector of labels.
##' @param par Vector of coefficients of split.
##' @return Gini impurity.
##' @author Johannes Tillil
gini_impurity <- function(dat, label, par) {
  select_idx <- as.matrix(dat) %*% par[2:length(par)] <= par[1]
  N1 <- sum(select_idx)
  N2 <- sum(!select_idx)
  if (N1 != 0 & N2 != 0) {
    gini <- 1-(
      N1/(N1+N2)*(
        (sum(label[select_idx] == "1")/N1)^2+
        (sum(label[select_idx] == "0")/N1)^2)+
      N2/(N1+N2)*(
        (sum(label[!select_idx] == "1")/N2)^2+
        (sum(label[!select_idx] == "0")/N2)^2)
    )
  } else if (N1 == 0) {
    gini <- 1-(
      (sum(label[!select_idx] == "1")/N2)^2+
      (sum(label[!select_idx] == "0")/N2)^2)
  } else if (N2 == 0) {
    gini <- 1-(
      (sum(label[select_idx] == "1")/N1)^2+
      (sum(label[select_idx] == "0")/N1)^2)
  }
  return(gini)
}

##' Compute Gini impurity of candidate linear splits in CART algorithm.
##' @title Gini impurity of candidate linear splits in CART algorithm.
##' @param dat Data frame of predictors.
##' @param label Vector of labels.
##' @param par Vector of coefficients of split.
##' @return Gini impurities.
##' @author Johannes Tillil
gini_impurity_CART <- function(dat, label, pure_coefs, value, candidates, varID, gamma) {
  coefs <- matrix(rep(pure_coefs, length(candidates)), ncol = length(candidates), byrow = FALSE)
  coefs[varID,] <- coefs[varID,] - candidates
  value <- value + candidates * gamma
  select_idx <- as.matrix(dat) %*% coefs <= drop(value)
  N1 <- colsums(select_idx)
  N2 <- colsums(!select_idx)
  gini <- 1-(
    N1/(N1+N2)*(
      (sapply(1:length(candidates), function(j) {sum(label[select_idx[,j]] == "1")})/N1)^2+
      (sapply(1:length(candidates), function(j) {sum(label[select_idx[,j]] == "0")})/N1)^2)+
    N2/(N1+N2)*(
      (sapply(1:length(candidates), function(j) {sum(label[!select_idx[,j]] == "1")})/N2)^2+
      (sapply(1:length(candidates), function(j) {sum(label[!select_idx[,j]] == "1")})/N2)^2)
  )
  sapply(
    1:length(candidates),
    function(candidateID) {
      if (N1[candidateID] == 0) {
        gini[candidateID] <<- 1-(
          (sum(label[!select_idx[,candidateID]] == "1")/N2[candidateID])^2+
          (sum(label[!select_idx[,candidateID]] == "0")/N2[candidateID])^2
        )
      }
      if (N2[candidateID] == 0) {
        gini[candidateID] <<- 1-(
          (sum(label[select_idx[,candidateID]] == "1")/N1[candidateID])^2+
          (sum(label[select_idx[,candidateID]] == "0")/N1[candidateID])^2
        )
      }
    }
  )
  return(gini)
}

##' Compute Batchwise gini impurity of candidate linear splits.
##' @title Batchwise gini impurity of candidate linear splits.
##' @param dat Data frame of predictors.
##' @param label Vector of labels.
##' @param par Vector of coefficients of split.
##' @return Gini impurities.
##' @author Johannes Tillil
gini_impurity_batch <- function(dat, label, candidates, varID) {
  coefs <- matrix(rep(numeric(ncol(dat)), length(candidates)), ncol = length(candidates), byrow = FALSE)
  coefs[varID,] <- 1
  value <- candidates
  select_idx <- as.matrix(dat) %*% coefs <= drop(value)
  N1 <- colsums(select_idx)
  N2 <- colsums(!select_idx)
  gini <- 1-(
    N1/(N1+N2)*(
      (sapply(1:length(candidates), function(j) {sum(label[select_idx[,j]] == "1")})/N1)^2+
      (sapply(1:length(candidates), function(j) {sum(label[select_idx[,j]] == "0")})/N1)^2)+
    N2/(N1+N2)*(
      (sapply(1:length(candidates), function(j) {sum(label[!select_idx[,j]] == "1")})/N2)^2+
      (sapply(1:length(candidates), function(j) {sum(label[!select_idx[,j]] == "1")})/N2)^2)
  )
  sapply(
    1:length(candidates),
    function(candidateID) {
      if (N1[candidateID] == 0) {
        gini[candidateID] <<- 1-(
          (sum(label[!select_idx[,candidateID]] == "1")/N2[candidateID])^2+
          (sum(label[!select_idx[,candidateID]] == "0")/N2[candidateID])^2
        )
      }
      if (N2[candidateID] == 0) {
        gini[candidateID] <<- 1-(
          (sum(label[select_idx[,candidateID]] == "1")/N1[candidateID])^2+
          (sum(label[select_idx[,candidateID]] == "0")/N1[candidateID])^2
        )
      }
    }
  )
  return(gini)
}

# ##' Compute Gini loss of linear split for Keras backend.
# ##' @title Gini loss of linear split
# ##' @param y_true Tensorflow tensor of true labels.
# ##' @param y_pred Tensorflow tensor of predicted labels.
# ##' @return Gini loss.
# ##' @author Johannes Tillil
# gini_loss <- function(y_true, y_pred) {
#   K <- backend()
#   y_pred <- (K$sign(y_pred-0.5)+1)/2
#   y_true <- K$transpose(y_true)-1
#   N1 <- K$sum(y_pred)$numpy()
#   N2 <- K$sum(1-y_pred)$numpy()
#   if (N1 != 0 & N2 != 0) {
#     1-(
#       N1/(N1+N2)*(
#         (K$dot(y_true, y_pred)/N1)^2+
#           (K$dot(y_true, 1-y_pred)/N1)^2) +
#         N2/(N1+N2)*(
#           (K$dot(1-y_true, y_pred)/N2)^2+
#             (K$dot(1-y_true, 1-y_pred)/N2)^2)
#     )
#   } else if (N1 == 0) {
#     1-(
#       (K$dot(1-y_true, y_pred)/N2)^2+
#         (K$dot(1-y_true, 1-y_pred)/N2)^2
#     )
#   } else if (N2 == 0) {
#     1-(
#       (K$dot(y_true, y_pred)/N1)^2+
#         (K$dot(y_true, 1-y_pred)/N1)^2
#     )
#   }
# }

go_to_parent <- function(child_nodeIDs, nodeID) {
  for (i in 1:length(child_nodeIDs)) {
    if (!is.null(child_nodeIDs[[i]])) {
      if (child_nodeIDs[[i]][1] == nodeID | child_nodeIDs[[i]][2] == nodeID) {
        break
      }
    }
  }
  return(i)
}
