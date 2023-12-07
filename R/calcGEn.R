#' @title calcGEn
#'
#' @description This function reads a CSV file containing single-cell gene expression data, performs Energy calculation using a Python script,
#' and returns the gene local network energy (GLNE) matrix and normalized scEnergy as a data frame.
#'
#' @param file_path Path to the input CSV file.
#' @param verbose Logical. If TRUE, print progress messages. Default is TRUE.
#' 
#' @return A list containing the gene local network energy matrix (GLNE) and a data frame with scEnergy and normalized scEnergy (scEn_Nor).
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- calcGEn("path/to/your/file.csv")
#' }
#'
#' @importFrom reticulate py_run_string
#'
#' @export


calcGEn <- function(file_path, verbose = TRUE) {
 
  # Validate if the file path exists
  if (!file.exists(file_path)) {
    stop("Error reading CSV file: ", file_path)
  }
  
  if (verbose) message(paste0(Sys.time(), ": Starting calculation."))

  py_code <- sprintf("

import numpy as np
import pandas as pd  

try:
    # Read CSV file data
    df = pd.read_csv('%s', index_col=0)
except Exception as e:
    raise RuntimeError(f'Invalid CSV file: {e}')

def calcGLNE(in_data):

    # Binarize data
    in_data = (in_data > 0).astype(int)

    n = in_data.shape[1]
    En = np.zeros((in_data.shape[0], n))

    for i in range(n):
        w = np.outer(in_data.iloc[:, i], in_data.iloc[:, i])
        np.fill_diagonal(w, 0)
        E1 = -1/2 * np.dot(in_data.iloc[:, i].T, w) * in_data.iloc[:, i]
        En[:, i] = E1

    En = pd.DataFrame(En, columns=in_data.columns, index=in_data.index)

    return En

# Call the Python function to calculate GLNE
GLNE_re = calcGLNE(df)
", file_path)

  py_run_string(py_code)
    
  GLNE <- py$GLNE_re
       
  scEnergy <- colSums(GLNE)
    
  scEn_Nor <- (scEnergy - min(scEnergy)) / (max(scEnergy) - min(scEnergy))
    
  scEn <- cbind(scEnergy,scEn_Nor)  
  scEn <- as.data.frame(scEn)
  colnames(scEn) <- c("scEnergy","scEn_Nor")
  
  if (verbose) message(paste0(Sys.time(), ": Complete GLNE and cell Energy calculation."))

  # Return Results
  return(list(GLNE = GLNE, scEn=scEn))   
}