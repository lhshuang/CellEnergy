#' @title find_terminal_states
#'
#' @description This function identifies terminal states in a numeric vector of normalized energy values for all cells.
#' The criteria for terminal state identification can be customized using various parameters.
#'
#' @param Nor_en A vector of normalized energy values for all cells. The name of the element is the index of the corresponding cell.
#' @param theta Threshold for ranking to identify terminal states. Default 0.85.
#' @param decreasing Logical value(TRUE or FALSE) indicating whether to rank in decreasing order, signifying a process in which cellular plasticity is reduced.
#' @param diffvar Logical value, parameter passed to PotencyFinder function.
#' @param maxPS Parameter passed to PotencyFinder function.
#'
#' @return A list containing cells with terminal states (`terminal_cells`) and filtered cells based on ranking (`candidate_cells`).
#'
#' @details
#' The function calculates ranks based on the specified criteria and identifies terminal states accordingly.
#' It utilizes the PotencyFinder function to further analyze and filter potential states.
#'
#' @examples
#' # Example usage:
#' find_terminal_states(Nor_en_vector)
#'
#' @importFrom stats qnorm
#'
#' @export

find_terminal_states <- function(Nor_en,theta=0.85,decreasing=TRUE,diffvar=TRUE,maxPS=3) {
  
   # Print a message indicating the identification of terminal states
   message(paste0(Sys.time(), ": Identification of terminal states."))
   
   # Check if the input is NULL or lacks named elements 
   if (is.null(Nor_en) || length(Nor_en) == 0 || is.null(names(Nor_en))) {
        stop("Error: 'Nor_en' is NULL or does not have named elements.")
    }   
    
   # Check if the PotencyFinder exists
   if (!exists("PotencyFinder", mode = "function")) {
     stop("Error: PotencyFinder function not found.")
   }
  
   # Check whether decreasing is logical value 
   if (!is.logical(decreasing) || length(decreasing) != 1) {
     stop("'decreasing' must be a single logical value (TRUE or FALSE).")
   }
    
  # Identify terminal states based on the specified criteria
  if (decreasing) {
      ranks <- 1-abs(Nor_en)
      names(ranks) <- names(Nor_en)
  
      # Calculate the cutoff point based on ranks
      cutoff <- qnorm(theta, mean=median(ranks), 
                  sd=median(abs(ranks - median(ranks))))
      
      # Select candidate terminal state cells with ranks above the cutoff point
      cell1 <- names(ranks)[ranks > cutoff] 
      
      # Identification of potential states
      En_filtered <- Nor_en[cell1]
      result <- suppressMessages(PotencyFinder(En_filtered,diffvar=diffvar,maxPS=maxPS))

      if (is.null(result$proM)) {
      stop("Error: Unexpected result from PotencyFinder.")
    }
      
      probS <- result$proM
      
      # Select cells with the lowest potential state
      colS <- ncol(probS)-1
      potS <- levels(factor(probS[, colS]))[1]   
      cell2 <- probS[probS[, colS] == potS,]
      cells <- rownames(cell2)
      
  } else {
      
      ranks <- abs(Nor_en)
      names(ranks) <- names(Nor_en)
  
 
      cutoff <- qnorm(theta, mean=median(ranks), 
                  sd=median(abs(ranks - median(ranks))))
      cell1 <- names(ranks)[ranks > cutoff] 
      
      En_filtered <-  Nor_en[cell1]
      result <- suppressMessages(PotencyFinder(En_filtered,diffvar=diffvar,maxPS=maxPS))
      
      if (is.null(result$proM)) {
      stop("Error: Unexpected result from PotencyFinder.")
    }
      
      probS <- result$proM
      
      # Select cells with the highest potential state
      colS <- ncol(probS) - 1
      potS <- tail(levels(factor(probS[, colS])), 1) 
      cell2 <- probS[probS[, colS] == potS,]
      cells <- rownames(cell2) 
  }
     message(paste0(Sys.time(), ": Done."))
    
     # Return a list containing cells with terminal states and candidate terminal cells based on ranking
     return(list(terminal_cells = cells, candidate_cells = cell1))
}