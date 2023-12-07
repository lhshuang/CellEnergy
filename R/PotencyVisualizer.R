#' @title PotencyVisualizer
#'
#' @description This function visualizes different potency states in 3D.
#'
#' @param data A list containing elements 'coordinates' and 'potencyState'.
#' @param numGrid Number of grid points for density calculations.
#' @param phi Angle specifying the colatitude of the plot. See more details in ?plot3D::persp3D
#' @param theta Angle specifying theazimuthal direction of the plot.
#' @param lphi Lighting angle for perspective plots.
#' @param k Number of colors to be used in the color palette.
#' @param colPersp Color palette for perspective plots.Suggest the minimum number of colors to meet 2*Number of potential states
#' @param colImage Color palette for image plots. Suggest using continuous colors
#' @param keyPersp List of parameters for perspective plot color key.
#' @param keyImage List of parameters for image plot color key.
#' @param width Width of the output plot.
#' @param height Height of the output plot.
#' @param plotPath Path to save the output plot. if plotPath=NULL; Default is the current working directory
#'
#' @return None.
#'
#' @examples
#' # Example usage:
#' PotencyVisualizer(myData)
#'
#' @import grDevices
#' @import graphics
#' @import plot3D
#' @import marray
#' @import MASS
#' @import marray
#'
#' @export

PotencyVisualizer <- function(data, numGrid = 50, phi = 20, theta = 20, lphi = 30, k = 15,
                             colPersp = NULL, colImage = NULL,
                             keyPersp = list(length = 0.2, width = 0.4, shift = 0.15, cex.axis = 0.7, cex.clab = 0.75),
                             keyImage = list(length = 0.2, width = 0.4, shift = -0.15, cex.axis = 0.7, cex.clab = 0.75),
                             width = 9, height = 3, plotPath = NULL) {
    
  # Check if the data structure meets the requirements
  if (!is.list(data) || !all(c("coordinates", "potencyState") %in% names(data))) {
    stop("Please provide a list containing elements'coordinates' and 'potencyState'!")
  }
       
  # Obtain the coordinates of cells
  coordinates <- data$coordinates
  component1 <- coordinates[, 1]
  component2 <- coordinates[, 2]

  # Calculate the density of all cells
  cellDensity <- kde2d(x = component1, y = component2, n = numGrid)

  potS <- data$potencyState
  idxP <- length(unique(potS))

  numRow <- ceiling(idxP / 3)
  numCol <- min(3, idxP)

  # Default color settings
  if (is.null(colPersp)) {
    colPersp <- colorRampPalette(c("red", "orange", "yellow2", "forestgreen", "blue"))(idxP * 2)
  }

  if (is.null(colImage)) {
    colImage <- maPalette(low = "white", mid = "#f09894", high = "red", k = k)
  }

  # Output the plot
  if (is.null(plotPath)) {
    plotPath <- getwd()
  }
  pdf(file.path(plotPath, "Potency_3D.pdf"), width = width, height = height)
  par(mar = c(2, 4, 2, 2))
  par(mfrow = c(numRow, numCol))

  idx1 <- 1
  idx2 <- 2

  for (i in 0:(idxP - 1)) {
    mainText <- if (i == 0) "Low potency state" else if (i == (idxP - 1)) "High potency state" else "Intermediate potency state"
    clabText <- if (i == 0) "LPS" else if (i == (idxP - 1)) "HPS" else "IPS"

    potency <- kde2d(x = component1[potS == (i + 1)],
                           y = component2[potS == (i + 1)],
                           n = 2 * numGrid)

    potency$z <- potency$z / max(potency$z)

    xlim <- c(min(range(cellDensity$x), range(potency$x)), max(range(cellDensity$x), range(potency$x)))
    ylim <- c(min(range(cellDensity$y), range(potency$y)), max(range(cellDensity$y), range(potency$y)))

    persp3D(x = potency$x, y = potency$y, z = (potency$z + (max(potency$z) * 0.5)) - (max(potency$z) * (0.2 * i)),
            main = mainText, xlim = xlim, ylim = ylim, zlim = c(-max(potency$z), max(potency$z)),
            phi = phi, theta = theta, col = maPalette(low = "lightgray", mid = colPersp[idx2], high = colPersp[idx1], k = k),
            colkey = keyPersp, lphi = lphi, clab = clabText, bty = "f", xlab = "Component 1", ylab = "Component 2")

    image3D(x = cellDensity$x, y = cellDensity$y, z = -max(potency$z), xlim = xlim, ylim = ylim,
            phi = phi, theta = theta, colvar = cellDensity$z, col = colImage, colkey = keyImage, clab = "Density(All)", add = TRUE)

    idx1 <- idx1 + 2
    idx2 <- idx2 + 2
  }
  dev.off()
  message(paste0(Sys.time(), ": Plotting completed."))
}
