#' @title findPath
#'
#' @description This function assigns fate paths to cells.
#' It calculates the category to which each cell belongs and creates a fate path.
#'
#' @param CSI A matrix containing cell state inference data. Output result of calcCSI function
#' @param En_Nor A data frame containing En_Nor data with columns 'scEn_Nor' and 'Class'.
#' @param decreasing A logical value indicating whether cellular plasticity should be decreasing (default is TRUE).
#'
#' @return A data frame On the basis of En_Nor, two columns of 'Path' and 'Pseudotime' are added, 
#' which correspond to the fate trajectory of the cell and the corresponding pseudotime value, respectively.
#'
#' @details
#' This function performs fate assignment by finding the category to which each cell
#' belongs based on the maximum values in the CSI matrix. It then creates a fate path
#' and calculates pseudotime for each cell.
#'
#' @examples
#' # Example usage:
#' result <- findPath(CSI_data, En_Nor_data)
#'
#' @seealso
#' \code{\link{check_columns}}
#'
#' @export

findPath <- function(CSI, En_Nor, decreasing = TRUE) {

  message(paste0(Sys.time(), ": Starting fate assignment."))
  
  # Check whether 'decreasing' is a logical value  
  if (!is.logical(decreasing) || length(decreasing) != 1) {
   stop("'decreasing' must be a single logical value (TRUE or FALSE).")
  }
 
  # Check that the En_Nor input is as expected  
  check_columns <- function(data, column_names) {
    missing_columns <- setdiff(column_names, colnames(data))
    if (length(missing_columns) > 0) {
      stop(paste("Data is missing columns:", paste(missing_columns, collapse = ", ")))
   }
  }
    
  check_columns(En_Nor, c("scEn_Nor", "Class"))

  if (!"Other" %in% En_Nor$Class) {
    stop("The 'Other' category is not found in the 'Class' column of 'En_Nor' data.")
  }
  
  # Check if row names of 'CSI' match the row names of 'En_Nor'
  if (!all(rownames(CSI) %in% rownames(En_Nor))) {
    stop("Error: Row names of 'CSI' must match the row names of 'En_Nor'.")
  }
  
  # Find the category to which each cell belongs
  column_categories <- apply(CSI, 1, function(x) {
    max_value <- max(x)
    colnames(CSI)[x == max_value]
  })

  Path <- list()

  # Get cell category
  categories <- unique(unlist(column_categories))

  # Traverse each category to find the corresponding cell
  for (category in categories) {
    values <- names(column_categories)[sapply(column_categories, function(x) category %in% x)]
    Path[[category]] <- if (length(values) > 0) values else NULL
  }

  # Filter data based on Class and prepare initial result
  result <- En_Nor[En_Nor$Class != "Other", ]
  result$Path <- result$Class
  result$Pseudotime <- rep(1, nrow(result))
                                              
  # Create fate path
  result_list <- list(result)
  for (i in seq_along(Path)) {
    group <- En_Nor[rownames(En_Nor) %in% Path[[i]], ]
    group$Path <- rep(names(Path[i]), nrow(group))
    group$Pseudotime <- if (decreasing) 1 - group$scEn_Nor else group$scEn_Nor
    result_list[[i + 1]] <- group
  }

  result <- bind_rows(result_list)
                                              
  message(paste0(Sys.time(), ": Done."))
  return(result)
}

