#' Covariance matrix visualisation
#'
#' @description Heatmap visualisation of the covariance matrix
#'    estimated within PPCA.
#'
#' @param ppcaOutput \code{list} -- the output object from running any
#'   of the PPCA functions in this package.
#'
#' @export
#'
#' @examples
ppca2Covplot <- function(ppcaOutput){
  paletteLength <- 50
  myColor <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                              "PuOr")))(paletteLength)
  myMax <- max(max(ppcaOutput$Sigma), abs(min(ppcaOutput$Sigma)))
  myBreaks <- c(seq(-myMax, 0, length.out=ceiling(paletteLength/2) + 1),
                seq(myMax/paletteLength, max(ppcaOutput$Sigma),
                    length.out=floor(paletteLength/2)))
  pheatmap::pheatmap(ppcaOutput$Sigma, color=myColor, breaks=myBreaks, fontsize =
                       8, cluster_rows = F, cluster_cols = F, legend = T, show_rownames = F,
                     show_colnames = F, border_color = NA)
}