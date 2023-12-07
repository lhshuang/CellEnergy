#' @title PotencyFinder
#' 
#' @description This function infers the discrete potency states of single cells and 
#' its distribution across the single cell population.
#' 
#' @param Nor_en
#' A vector of normalized energy values for all cells. The name of the element is the index of the corresponding cell.
#' 
#' @param pheno
#' (Optional) A vector representing phenotypic information for each cell.
#' 
#' @param diffvar
#' A logical value (TRUE or FALSE) indicating whether to assume different variances for clusters in the Gaussian Mixture Model.Default is TRUE.
#' 
#' @param maxPS
#' An integer specifying the maximum number of potency states to infer.
#' 
#' @param verbose
#' A logical value (TRUE or FALSE) indicating whether to display detailed output information.
#' 
#' @param path
#' (Optional) A character string specifying the path to save output plots and model information.
#' 
#' @return 
#' A list containing the Bayesian Information Criterion (BIC), probability matrix (proM), 
#' inferred potency states (potS), distribution of potency states across phenotypes (distPS), 
#' probability matrix for each phenotype (probPS), and normalised Heterogeneous Index for each phenotype (HI).
#' 
#' @examples
#' # Example usage:
#' result <- PotencyFinder(Nor_en_vector, pheno_vector, diffvar = TRUE, maxPS = 5, verbose = TRUE, path = "output_path")
#' 
#' @import mclust
#'
#' @export

PotencyFinder <- function(Nor_en, pheno = NULL, diffvar = TRUE,maxPS = 5,verbose=FALSE,path=NULL){
  
  
    # fit Gaussian Mixture Model for potency inference
    message(paste0(Sys.time(), ": Fitting GMM to logit-transformed Energy values."))
    
    if (is.null(Nor_en) || length(Nor_en) == 0 || is.null(names(Nor_en))) {
        stop("Error: 'Nor_en' is NULL or does not have named elements.")
    }
     
    logitEn <- log2(1+Nor_en)
    
    # Check whether diffvar is logical value 
    if (!is.logical(diffvar) || length(diffvar) != 1) {
        stop("'diffvar' must be a single logical value (TRUE or FALSE).")
    }
    
    # Check whether verbose is logical value 
    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("'verbose' must be a single logical value (TRUE or FALSE).")
    }
    
    if(diffvar){ 
     # default assumes different variance for clusters
        GMM_en <- Mclust(logitEn, G = seq_len(maxPS))  # G个clusters variance different
    } else {
        GMM_en <- Mclust(logitEn, G = seq_len(maxPS), modelNames = c("E"))
    }
    
    
    # Output plots
    if (!is.null(path)) {
        pdf(file.path(path, "plot.pdf"))
        plot(GMM_en, what = "BIC", color = "red")
        plot(GMM_en, what = "classification")
        plot(GMM_en, what = "uncertainty")
        dev.off()
        print(paste0(Sys.time(), ": Plots generated."))
    }
    
    # Output Model Information
    if (!is.null(path) && verbose) {
        parameters_output <- capture.output(GMM_en$parameters)
        summary_output <- capture.output(summary(GMM_en))
        cat(parameters_output, summary_output, file = file.path(path, "GMM_summary.txt"), sep = "\n")
        print(paste0(Sys.time(), ": Model Output."))
    }
        
    
    bic <- GMM_en$BIC
    probability <- GMM_en$z

    potS <- GMM_en$classification
    nPS <- length(unique(potS))
    print(paste0(Sys.time(), ": Identified ",nPS," potency states."))
    
    
    ord <- as.numeric(levels(factor(potS)))
    ordpotS <- match(potS,ord)
    
    GMM_en$classification <- ordpotS   
   
    for (i in seq_len(nPS)) {
        names(ordpotS)[ordpotS == i] <- paste("PS", i, sep = "")
    }
    
    probability <- cbind(probability,GMM_en$classification,names(ordpotS))
    rownames(probability) <- names(Nor_en)
    
    # calculate the heterogeneity index for each phenotype category
    if(!is.null(pheno)){
      nPH <- length(unique(pheno))
      distPSph <- table(pheno,ordpotS) 
      message(paste0(Sys.time(), ": Calculate the Heterogeneous index for each phenotype category."))
      probPSph <- distPSph / rowSums(distPSph)
  
      het <- sapply(seq_len(nPH), function(ph) {
          prob <- probPSph[ph, ]
          idx <- prob > 0
         -sum(prob[idx] * log(prob[idx])) / log(nPS)   
      })
      
        names(het) <- rownames(probPSph)
        
        message(paste0(Sys.time(), ": Complete."))
    } else {
        distPSph=NULL
        probPSph=NULL
        het=NULL
    }
    
    return(list(BIC=bic, proM =probability, potS = ordpotS, distPS=distPSph, probPS=probPSph, HI=het))
}