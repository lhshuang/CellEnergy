#' @title calcCSI
#'
#' @description This function calculates the Cell Similarity Index (CSI) based on the input terminal data and gene local network energy (GLNE) data.
#'
#' @param ter_data A character vector with named elements representing terminal cell types and remaining cells (Other).
#' @param GLNE_data Gene local network energy matrix where columns represent cells and rows represent genes.
#' @param normalize A logical value indicating  whether the GLNE_data needs to be processed  (default is TRUE).
#'
#' @return A matrix representing the Cell Similarity Index (CSI) between remaining cells and terminal cell types.
#' Columns represent terminal cell types and rows represent remaining cells.
#'
#' @examples
#' # Example usage:
#' ter_data <- c("Cell1" = "Type1", "Cell2" = "Type2", "Cell3" = "Other") 
#' GLNE_data <- matrix(rnorm(100), nrow = 10, ncol = 10)  # Example GLNE data
#' result <- calcCSI(ter_data, GLNE_data, normalize =TRUE) # Names in ter_data must match the column names of GLNE_data.
#'
#' @import philentropy
#'
#' @export

calcCSI <- function(ter_data, GLNE_data,normalize=TRUE) {
    
  message(paste0(Sys.time(), ": Starting CSI calculation."))
    
  # Check if the input is factor or lacks named elements 
  if (!is.character(ter_data) || is.null(names(ter_data))) {
     stop("Error: 'ter_data' must be a character vector and have named elements.")
   } 
    
  # Check if names in 'ter_data' match the column names of 'GLNE_data'
  if (!all(names(ter_data) %in% colnames(GLNE_data))) {
     stop("Error: Names in 'ter_data' must match the column names of 'GLNE_data'.")
   }

  # Check whether 'normalize' is a logical value  
  if (!is.logical(normalize) || length(normalize) != 1) {
     stop("'normalize' must be a single logical value (TRUE or FALSE).")
   } 
    
  # Get levels of the categorical variable
  terC <- levels(factor(ter_data))
    
  # Preprocess input data  
  if (normalize) {
    GLNE0 <- apply(GLNE_data, 2, function(x) {
      if (sum(x) == 0) {
        return(log(x / sum(x - 0.001) * 1e4 + 1))
      } else {
        return(log(x / sum(x) * 1e4 + 1))
      }  
    })
    tGLNE <- apply(GLNE0, 1, function(x) (x - mean(x)) / sd(x))
    GLNE_data1 <- t(tGLNE)
   } else {
    tGLNE <- t(GLNE_data) 
    GLNE_data1 <- GLNE_data
 }
    
  # Create subsets and name them based on levels(terminal state cell type)
  subsets <- setNames(lapply(terC, function(level) ter_data[ter_data == level]), terC)

  # Extract terminal state cells
  filtered_terC <- terC[terC != 'Other']
  
  if (length(terC) > 2){
      # Calculate different terminal states
      terSM <- lapply(filtered_terC, function(level) {
        subset_data <- tGLNE[names(subsets[[level]]), ]
        colMeans(subset_data)
      })
      
      terSM <- lapply(filtered_terC, function(level) colMeans(tGLNE[names(subsets[[level]]), ]))

      # Extract the GLNE matrix corresponding to the remaining cells
      remC  <- GLNE_data1[, names(subsets[["Other"]])]

      # Initialize distance matrix and correlation matrix
      GD <- matrix(0, nrow = ncol(remC), ncol = length(terSM))
      PC <- matrix(0, nrow = ncol(remC), ncol = length(terSM))
      
      # Calculate the Gore distance between each remaining cell and the terminal state
      for (i in 1 : ncol(remC)) {
        for (j in 1 : length(terSM)) {
          GD[i, j] <- philentropy::distance(rbind(remC[, i], terSM[[j]]), method = "gower", mute.message = TRUE)
          }
      } 
      
      # Calculate the Pearson correlation between each remaining cell and the terminal state
      PC <- sapply(terSM, function(mean_vector) cor(remC, mean_vector, method = "pearson", use = "pairwise"))
                   
      # Set rowname and colname
      rownames(GD) <- colnames(remC)
      colnames(GD) <- filtered_terC
    
      rownames(PC) <- colnames(remC)
      colnames(PC) <- filtered_terC
                   
      # Calculate CSI 
      CSI_Ori <- (1 - GD) * PC            
      CSI_Nor <- t(apply(CSI_Ori, 1, function(row) (row - min(row)) / (max(row) - min(row))))
      CSI <- CSI_Nor/rowSums(CSI_Nor)
    } else {
      remC <- names(subsets[["Other"]])
      CSI <- matrix(1, nrow = length(remC), ncol = 1)
      row.names(CSI) <- remC
      colnames(CSI) <- filtered_terC
    }
    message(paste0(Sys.time(), ": Calculation completed."))
    return(CSI)
  }