##' Generates random networks of genes, their associations, expression data and
##' binary labels where the effect measures of individual genes can be set by
##' the user.
##' 
##' @title genGeneNetworkData
##' @param num_networks Integer, number of networks to generate.
##' @param num_genes Integer, number of genes per network.
##' @param num_modules Integer, number of modules per network. Can also be NULL for random number of modules.
##' @param num_observations Integer, number of expression data observations per network.
##' @param num_causal_modules Integer, how many modules should be causal for the phenotype.
##' @param prop_causal_genes Float, proportion of genes in each causal module that should be causal for the phenotype. If "causal_genes_randomly_distributed == TRUE" refers to proportion of all genes that should be causal.
##' @param effect_size Float, standardized effect size of causal genes in the network.
##' @param effect_intercept Float, standardized intercept effect of all genes in the network.
##' @param causal_genes_randomly_distributed Boolean, if TRUE causal genes will be sampled randomly from all available genes.
##' @param num_threads Integer, number of cores to parallelize on.
##' @param seed Integer, initial seed for the L'Ecuyer-CMRG random number streams.
##' @param effect_type Character, type of gene-effect on phenotype. Default is "linear". One of "linear", "quadratic".
##' @examples
##' \donttest{
##' library(simpleRFNetwork)
##' 
##' 
##' }
##' 
##' @author Johannes Tillil
##' @import stats
##' @import parallel
##' @export
genGeneNetworkData <- function(
  num_networks,
  num_genes,
  num_modules = NULL,
  num_observations,
  num_causal_modules,
  prop_causal_genes,
  effect_size = 1,
  effect_intercept = -1,
  effect_type = "linear",
  effect_error_sd = 0,
  causal_genes_randomly_distributed = FALSE,
  num_threads = 1,
  seed = 1
) {

  ## Set up parallel reproducibility
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  mc.reset.stream()
  
  ## Start parallel computing
  return(mclapply(
    1:num_networks,
    function(i) {
      ## Generate networks, modules and associations between genes
      if (is.null(num_modules)) {
        rn <- random_network(num_genes)
      } else {
        rn <- random_network(num_genes, num_modules)
      }
      rn <- gen_partial_correlations(rn)
      exprdat <- scale(log(gen_rnaseq(num_observations, rn)$x + 1))
      modules <- lapply(rn$modules, function(x){x$nodes})
      num_modules <- length(rn$modules)

      ## Sample causal modules and genes and set gene effects accordingly
      if (causal_genes_randomly_distributed) {
        ## Sample causal genes
        causal_genes <- sample(
          1:num_genes,
          ceiling(prop_causal_genes*num_genes),
          replace = FALSE
        )
        ## Count causal modules
        causal_modules <- NULL
        sapply(
          1:length(modules),
          function(j) {
            if (sum(causal_genes %in% modules[[j]]) > 0) {
              causal_modules <<- c(causal_modules, j)
            }
          }
        )
        ## Save sampled values
        causal_genes <- unique(causal_genes)
        effects <- numeric(num_genes)
        effects[causal_genes] <- effect_size
        effects <- effects
        causal_modules <- causal_modules
        causal_genes <- causal_genes
      } else if (num_causal_modules > 0 & prop_causal_genes > 0) {
        ## Sample causal modules
        causal_modules <- sample(
          x = 1:num_modules,
          size = min(num_causal_modules, num_modules),
          replace = FALSE)
        ## Sample causal genes
        causal_genes <- NULL
        sapply(
          1:min(num_causal_modules, num_modules),
          function(j) {
            ## Read required amount of genes
            num_required_genes <- ceiling(prop_causal_genes*length(modules[[causal_modules[j]]]))
            ## Sample first causal gene
            sampled_genes <- sample(
              x = 1:length(modules[[causal_modules[j]]]),
              size = 1,
              replace = FALSE
            )
            ## Read adjacency matrix for module
            adj_mat <- get_adjacency_matrix(rn)[modules[[causal_modules[j]]], modules[[causal_modules[j]]]]
            ## Search for associated genes in the module
            while (length(sampled_genes) < num_required_genes) {
              ## For all candidates not yet added to causal genes
              for (candidate_id in (1:nrow(adj_mat))[-sampled_genes]) {
                ## If still required AND connected to sampled_genes
                if (
                  length(sampled_genes) < num_required_genes &
                  sum(adj_mat[candidate_id, sampled_genes]) > 0
                ) {
                  sampled_genes <- c(sampled_genes, candidate_id)
                }
              }
            }
            causal_genes <<- c(causal_genes, modules[[causal_modules[j]]][sampled_genes])
          }
        )
        ## Save sampled values
        causal_genes <- unique(causal_genes)
        effects <- numeric(num_genes)
        effects[causal_genes] <- effect_size
        effects <- effects
        causal_modules <- causal_modules
        causal_genes <- causal_genes
      } else {
        causal_modules <- NULL
        causal_genes <- NULL
        effects <- numeric(num_genes)
      }
      
      ## Sample phenotype from gene effects and combine phenotype and expression data into one data frame for training
      if (effect_error_sd > 0) {
        effects <- effects + rnorm(length(effects), 0, effect_error_sd)
      }
      if (effect_type == "linear") {
        probs <- 1/(1 + exp(-as.matrix(exprdat) %*% effects - effect_intercept))
      } else if (effect_type == "quadratic") {
        probs <- 1/(1 + exp(-(sign(as.matrix(exprdat) %*% effects) * (as.matrix(exprdat) %*% effects)^2) - effect_intercept))
      }  else if (effect_type == "cubic") {
        probs <- 1/(1 + exp(-((as.matrix(exprdat) %*% effects)^3) - effect_intercept))
      }
      res <- data.frame(pheno = as.factor(sapply(
        1:num_observations,
        function(i) {
          sample(
            x = c(0,1),
            size = 1,
            prob = c(1-probs[i], probs[i]))
        }
      )))
      res <- cbind(res, data.frame(exprdat))
      return(list(
        data = res, 
        modules = modules,
        causal_modules = causal_modules,
        causal_genes = causal_genes,
        effects = effects
      ))
    },
    mc.cores = num_threads
  ))
}
