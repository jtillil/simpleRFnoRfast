univariate_split <- function(data_values, response) {
  ## Initiate
  Gini_impurity_start <- 9999
  best_val <- 0
  best_varID <- 0
  
  ## Iterate over all variables
  sapply(1:ncol(data_values), function(varID) {
    ## Read value candidates
    val_candidates <- unique(data_values[,varID])
    sapply(val_candidates, function(val) {
      ## Compute new Gini impurity
      coefficients_start <- numeric(ncol(data_values))
      coefficients_start[varID] <- 1
      Gini_impurity_val <- gini_impurity(data_values,
                                         response,
                                         c(val, coefficients_start))
      
      ## Compare to current best Gini impurity and set value, varID if smaller
      if (Gini_impurity_val < Gini_impurity_start) {
        Gini_impurity_start <<- Gini_impurity_val
        best_val <<- val
        best_varID <<- varID
      }
    })
  })
  
  ## Coerce best uni-variate split into coefficients and value
  coefficients <- numeric(ncol(data_values))
  coefficients[best_varID] <- 1
  value <- best_val
  
  return(c(value, coefficients))
}

univariate_split_fast <- function(data_values, response) {
  ## Set fraction of subset variables
  nu <- 0.1
  minN <- 100
  
  ## Initiate
  Gini_impurity_start <- 9999
  best_val <- 0
  best_varID <- 0
  
  ## Iterate over all variables
  sapply(1:ncol(data_values), function(varID) {
    ## Sample value candidates
    unique_col <- unique(data_values[,varID])
    val_candidates <- sample(unique_col, min(max(round(nu*nrow(data_values)), minN), length(unique_col)))
    sapply(val_candidates, function(val) {
      ## Compute new Gini impurity
      coefficients_start <- numeric(ncol(data_values))
      coefficients_start[varID] <- 1
      Gini_impurity_val <- gini_impurity(data_values,
                                         response,
                                         c(val, coefficients_start))
      
      ## Compare to current best Gini impurity and set value, varID if smaller
      if (Gini_impurity_val < Gini_impurity_start) {
        Gini_impurity_start <<- Gini_impurity_val
        best_val <<- val
        best_varID <<- varID
      }
    })
  })
  
  ## Coerce best uni-variate split into coefficients and value
  coefficients <- numeric(ncol(data_values))
  coefficients[best_varID] <- 1
  value <- best_val
  
  return(c(value, coefficients))
}

SVM <- function(data_values, response) {
  ## Calculate SVM plane
  svmfit <- svm(y=response,
                x=data_values,
                kernel="linear",
                scale=FALSE)
  
  ## Read coefficients and value
  coefficients <- drop(t(svmfit$coefs)%*%as.matrix(data_values)[svmfit$index,])
  value <- svmfit$rho
  
  return(c(value, coefficients))
}

LDA <- function(data_values, response, mat) {
  ## Calculate class means
  mean0 <- colmeans(as.matrix(data_values[response == 0,]))
  mean1 <- colmeans(as.matrix(data_values[response == 1,]))
  
  ## Calculate coefficients and value
  coefficients <- spdinv(mat) %*% (mean1 - mean0)
  value <- sum(coefficients * (0.5*(mean1 + mean0)))
  
  return(c(value, coefficients))
}

gini_optim <- function(data_values, response) {
  ## Find first split as best uni-variate split
  ## Set fraction of subset variables
  nu <- 0.1
  minN <- 100
  
  ## Initiate
  Gini_impurity_start <- 9999
  best_val <- 0
  best_varID <- 0
  
  ## Iterate over all variables
  sapply(1:ncol(data_values), function(varID) {
    ## Sample value candidates
    unique_col <- unique(data_values[,varID])
    val_candidates <- sample(unique_col, min(max(round(nu*nrow(data_values)), minN), length(unique_col)))
    sapply(val_candidates, function(val) {
      ## Compute new Gini impurity
      coefficients_start <- numeric(ncol(data_values))
      coefficients_start[varID] <- 1
      Gini_impurity_val <- gini_impurity(data_values,
                                         response,
                                         c(val, coefficients_start))
      
      ## Compare to current best Gini impurity and set value, varID if smaller
      if (Gini_impurity_val < Gini_impurity_start) {
        Gini_impurity_start <<- Gini_impurity_val
        best_val <<- val
        best_varID <<- varID
      }
    })
  })
  
  ## Coerce best uni-variate split into parameters
  coefficients <- numeric(ncol(data_values))
  coefficients[best_varID] <- 1
  param <- c(best_val, coefficients)
  
  ## Calculate Gini-optimal plane with univariate initial split
  par <- optim(
    par = param,
    fn = function(par) {
      return(gini_impurity(data_values,
                           response,
                           par))
    },
    method="Nelder"
  )$par
  
  ## Read coefficients and value
  coefficients <- par[2:length(par)]
  value <- par[1]
  
  return(c(value, coefficients))
}

stoch_optim <- function(data_values, response) {
  ## Turn on eager execution
  tf$config$run_functions_eagerly(TRUE)
  
  ## Create, compile and fit Keras model
  model <- keras_model_sequential() %>%
    layer_dense(
      units = 1,
      activation = "sigmoid",
      input_shape = dim(data_values)[[2]],
      kernel_initializer='glorot_uniform',
      bias_initializer='glorot_uniform'
    )
  model %>%
    compile(
      optimizer = tf$keras$optimizers$Adam(learning_rate = 0.1),
      # optimizer = tf$keras$optimizers$SGD(learning_rate = 0.1),
      loss = "binary_crossentropy",
      run_eagerly = TRUE
    )
  model %>%
    fit(as.matrix(data_values),
        as.numeric(response)-1,
        epochs = 100,
        verbose = 0,
        callbacks = c(callback_early_stopping(
          monitor = "loss",
          min_delta = 0.01,
          mode = "min",
          patience = 10
        ))
    )
  
  ## Read coefficients and value 
  coefficients <- as.numeric(model$layers[[1]]$weights[[1]]$numpy())
  value <- -as.numeric(model$layers[[1]]$weights[[2]]$numpy())
  
  return(c(value, coefficients))
}

CART <- function(IQR_data_values, data_values, response) {
  ## Find first split as best uni-variate split
  ## Initiate
  Gini_impurity_start <- 9999
  best_val <- 0
  best_varID <- 0
  
  ## Iterate over all variables
  sapply(1:IQR_data_values$ncol, function(varID) {
    ## Read value candidates
    val_candidates <- unique(IQR_data_values$column(varID))
    sapply(val_candidates, function(val) {
      ## Compute new Gini impurity
      coefficients_start <- numeric(IQR_data_values$ncol)
      coefficients_start[varID] <- 1
      Gini_impurity_val <- gini_impurity(IQR_data_values$data,
                                         response,
                                         c(val, coefficients_start))
      
      ## Compare to current best Gini impurity and set value, varID if smaller
      if (Gini_impurity_val < Gini_impurity_start) {
        Gini_impurity_start <<- Gini_impurity_val
        best_val <<- val
        best_varID <<- varID
      }
    })
  })
  
  ## Coerce best uni-variate split into coefficients and value
  coefficients <- numeric(IQR_data_values$ncol)
  coefficients[best_varID] <- 1
  value <- best_val
  
  ## Compute starting Gini impurity
  Gini_impurity_nplus1 <- gini_impurity(IQR_data_values$data,
                                        response,
                                        c(value, coefficients))
  Gini_impurity_n <- 9999
  
  ## Set convergence threshold
  epsilon <- 0.01
  
  ## Perform updates until improvement below threshold
  while ((Gini_impurity_n - Gini_impurity_nplus1) > epsilon) {
    ## Set Gini impurity of last cycle
    Gini_impurity_n <- Gini_impurity_nplus1
    
    ## Cycle through all variables and search for an improved split by varying their coefficient
    sapply(1:IQR_data_values$ncol,
           function(varID) {
             ## Compute current split values for all observations
             v <- as.matrix(IQR_data_values$data) %*% coefficients
             
             ## For gamma equals -0.25, 0, 0.25
             sapply(c(-0.25, 0, 0.25), function(gamma) {
               ## Compute candidates
               u <- (v - value) / (IQR_data_values$column(varID) + gamma)
               
               ## For every candidate
               sapply(u, function(u_n) {
                 ## Convert to candidate coefficients and value
                 coefficients_u <- coefficients
                 coefficients_u[varID] <- coefficients_u[varID] - u_n
                 value_u <- value + u_n * gamma
                 
                 if (gini_impurity(IQR_data_values$data,
                                   response,
                                   c(value_u, coefficients_u)
                 ) < Gini_impurity_nplus1) {
                   
                   ## Update coefficients, value, split values
                   coefficients <<- coefficients_u
                   value <<- value_u
                   v <<- as.matrix(IQR_data_values$data) %*% coefficients
                   
                   ## Calculate new Gini impurity
                   Gini_impurity_nplus1 <<- gini_impurity(IQR_data_values$data,
                                                          response,
                                                          c(value, coefficients))
                 }
               })
             })
           }
    )
  }
  
  ## Read IQR scaling parameters
  IQR_vals <- c()
  sapply(1:ncol(data_values), function(varID) {
    IQR_vals <<- c(IQR_vals, IQR(data_values[,varID]))
  })
  Mean_vals <- colmeans(as.matrix(data_values))
  
  ## Rescale coefficients and value
  axis_points <- value/coefficients*IQR_vals
  translated_axis_points <- NULL
  for (varID in 1:length(coefficients)) {
    translated_axis_points[varID] <- axis_points[varID] + Mean_vals[varID] + sum(axis_points[varID]/axis_points[-varID] * Mean_vals[-varID])
  }
  coefficients <- 1/translated_axis_points
  value <- 1

  ## Set NaN coefficients to 0
  coefficients[is.nan(coefficients)] <- 0
  
  return(c(value, coefficients))
}

# CART_fast <- function(IQR_data_values, data_values, response) {
#   ## Find first split as best uni-variate split
#   ## Set fraction of subset variables
#   nu <- 0.1
  
#   ## Initiate
#   Gini_impurity_start <- 9999
#   best_val <- 0
#   best_varID <- 0
  
#   sapply(1:IQR_data_values$ncol, function(varID) {
#     ## Sample value candidates
#     unique_col <- unique(IQR_data_values$column(varID))
#     val_candidates <- sample(unique_col, round(nu*length(unique_col)))
#     sapply(val_candidates, function(val) {
#       ## Compute new Gini impurity
#       coefficients_start <- numeric(IQR_data_values$ncol)
#       coefficients_start[varID] <- 1
#       Gini_impurity_val <- gini_impurity(
#         IQR_data_values$data,
#         response,
#         c(val, coefficients_start)
#       )
      
#       ## Compare to current best Gini impurity and set value, varID if smaller
#       if (Gini_impurity_val < Gini_impurity_start) {
#         Gini_impurity_start <<- Gini_impurity_val
#         best_val <<- val
#         best_varID <<- varID
#       }
#     })
#   })
  
#   ## Coerce best uni-variate split into coefficients and value
#   coefficients <- numeric(IQR_data_values$ncol)
#   coefficients[best_varID] <- 1
#   value <- best_val
  
#   ## Compute starting Gini impurity
#   Gini_impurity_nplus1 <- gini_impurity(
#     IQR_data_values$data,
#     response,
#     c(value, coefficients))
#   Gini_impurity_n <- 9999
  
#   ## Set convergence threshold
#   epsilon <- 0.001
  
#   ## Perform updates until improvement below threshold
#   while ((Gini_impurity_n - Gini_impurity_nplus1) > epsilon) {
#     ## Set Gini impurity of last cycle
#     Gini_impurity_n <- Gini_impurity_nplus1
    
#     ## Cycle through all variables and search for an improved split by varying their coefficient
#     sapply(
#       1:IQR_data_values$ncol,
#       function(varID) {
#         ## Sample subset
#         subset_data_values <- IQR_data_values$data[sample(1:nrow(data_values), round(nu*nrow(data_values))),]

#         ## Compute current split values for subset observations
#         v <- as.matrix(subset_data_values) %*% coefficients
        
#         ## For gamma equals -0.25, 0, 0.25
#         sapply(c(-0.25, 0, 0.25), function(gamma) {
#           ## Compute candidates
#           u <- (v - value) / (subset_data_values[,varID] + gamma)
          
#           ## For every candidate
#           sapply(u, function(u_n) {
#             ## Convert to candidate coefficients and value
#             coefficients_u <- coefficients
#             coefficients_u[varID] <- coefficients_u[varID] - u_n
#             value_u <- value + u_n * gamma
            
#             if (gini_impurity(
#               IQR_data_values$data,
#               response,
#               c(value_u, coefficients_u)
#             ) < Gini_impurity_nplus1) {
              
#               ## Update coefficients, value, split values
#               coefficients <<- coefficients_u
#               value <<- value_u
#               v <<- as.matrix(IQR_data_values$data) %*% coefficients
              
#               ## Calculate new Gini impurity
#               Gini_impurity_nplus1 <<- gini_impurity(
#                 IQR_data_values$data,
#                 response,
#                 c(value, coefficients)
#               )
#             }
#           })
#         })
#       }
#     )
#   }
  
#   ## Read IQR scaling parameters
#   IQR_vals <- c()
#   sapply(1:ncol(data_values), function(varID) {
#     IQR_vals <<- c(IQR_vals, IQR(data_values[,varID]))
#   })
#   Mean_vals <- colmeans(as.matrix(data_values))
  
#   ## Rescale coefficients and value
#   axis_points <- value/coefficients*IQR_vals
#   translated_axis_points <- NULL
#   for (varID in 1:length(coefficients)) {
#     translated_axis_points[varID] <- axis_points[varID] + Mean_vals[varID] + sum(axis_points[varID]/axis_points[-varID] * Mean_vals[-varID])
#   }
#   coefficients <- 1/translated_axis_points
#   value <- 1

#   ## Set NaN coefficients to 0
#   coefficients[is.nan(coefficients)] <- 0
  
#   return(c(value, coefficients))
# }

CART_fast <- function(IQR_data_values, data_values, response) {
  ## Find first split as best uni-variate split
  ## Set fraction of subset variables
  nu <- 0.1
  minN <- 100
  
  ## Initiate
  Gini_impurity_start <- 9999
  best_val <- 0
  best_varID <- 0
  
  sapply(1:IQR_data_values$ncol, function(varID) {
    ## Sample value candidates
    unique_col <- unique(IQR_data_values$column(varID))
    val_candidates <- sample(unique_col, min(max(round(nu*nrow(data_values)), minN), length(unique_col)))
    sapply(val_candidates, function(val) {
      ## Compute new Gini impurity
      coefficients_start <- numeric(IQR_data_values$ncol)
      coefficients_start[varID] <- 1
      Gini_impurity_val <- gini_impurity(
        IQR_data_values$data,
        response,
        c(val, coefficients_start)
      )
      
      ## Compare to current best Gini impurity and set value, varID if smaller
      if (Gini_impurity_val < Gini_impurity_start) {
        Gini_impurity_start <<- Gini_impurity_val
        best_val <<- val
        best_varID <<- varID
      }
    })
  })
  
  ## Coerce best uni-variate split into coefficients and value
  coefficients <- numeric(IQR_data_values$ncol)
  coefficients[best_varID] <- 1
  value <- best_val
  
  ## Compute starting Gini impurity
  Gini_impurity_nplus1 <- gini_impurity(
    IQR_data_values$data,
    response,
    c(value, coefficients)
  )
  Gini_impurity_n <- 9999
  
  ## Set convergence threshold
  epsilon <- 0.001
  
  ## Perform updates until improvement below threshold
  while ((Gini_impurity_n - Gini_impurity_nplus1) > epsilon) {
    ## Set Gini impurity of last cycle
    Gini_impurity_n <- Gini_impurity_nplus1

    ## Sample subset
    subset_data_values <- IQR_data_values$data[sample(1:nrow(data_values), min(max(round(nu*nrow(data_values)), minN), nrow(data_values))),]

    ## Cycle through all variables and search for an improved split by varying their coefficient
    sapply(
      1:IQR_data_values$ncol,
      function(varID) {
        ## Compute current split values for subset observations
        v <- as.matrix(subset_data_values) %*% coefficients

        ## For gamma equals -0.25, 0, 0.25
        best_u <- sapply(c(-0.25, 0, 0.25), function(gamma) {
          ## Compute candidates
          u <- (v - value) / (subset_data_values[,varID] + gamma)
          u[is.nan(u)] <- 0
          
          ## Calculate gini impurities for candidates
          impurities <- gini_impurity_CART(
            dat = IQR_data_values$data,
            label = response,
            pure_coefs = coefficients,
            value = value,
            candidates = u,
            varID = varID,
            gamma = gamma
          )

          ## Extract best candidate
          best_idx <- which.min.random(impurities)
          c(impurities[best_idx], u[best_idx])
        })

        ## Extract best candidate
        best_gamma <- c(-0.25, 0, 0.25)[which.min.random(best_u[1,])]
        best_impurity <- best_u[1, which.min.random(best_u[1,])]
        best_u <- best_u[2, which.min.random(best_u[1,])]

        if (best_impurity < Gini_impurity_nplus1) {
          ## Update coefficients, value, split values, gini_impurity
          coefficients[varID] <<- coefficients[varID] - best_u
          value <<- value + best_u * best_gamma
          new_coefficients <<- TRUE
          Gini_impurity_nplus1 <<- best_impurity
        }
      }
    )
  }

  ## Read IQR scaling parameters
  IQR_vals <- c()
  sapply(1:ncol(data_values), function(varID) {
    IQR_vals <<- c(IQR_vals, IQR(data_values[,varID]))
  })
  Mean_vals <- colmeans(as.matrix(data_values))
  
  ## Rescale coefficients and value
  axis_points <- value/coefficients*IQR_vals
  translated_axis_points <- NULL
  for (varID in 1:length(coefficients)) {
    translated_axis_points[varID] <- axis_points[varID] + Mean_vals[varID] + sum(axis_points[varID]/axis_points[-varID] * Mean_vals[-varID])
  }
  coefficients <- 1/translated_axis_points
  value <- 1

  ## Set NaN coefficients to 0
  coefficients[is.nan(coefficients)] <- 0
  
  return(c(value, coefficients))
}
