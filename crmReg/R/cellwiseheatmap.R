cellwiseheatmap <- function (cellwiseoutliers, data,
                             col = c("blue", "lightgray", "red"),
                             notecol.outlier = "white", notecol.clean = "black", notecex = 1,
                             margins = c(9.5, 14), lhei = c(0.5, 15), lwid = c(0.1, 3.5),
                             sepcolor = "white", sepwidth = c(0.01, 0.01)) {
  n <- nrow(cellwiseoutliers)
  p <- ncol(cellwiseoutliers)
  rownames(cellwiseoutliers) <- rownames(data)
  colnames(cellwiseoutliers) <- colnames(data)
  notecol <- ifelse(c(t(cellwiseoutliers[n:1, ])) != 0, notecol.outlier, notecol.clean)
  heatmap.2(as.matrix(cellwiseoutliers),
            scale        = "none",
            col          = col,
            cellnote     = data,
            notecol      = notecol,
            notecex      = notecex,
            Rowv         = FALSE,
            Colv         = FALSE,
            dendrogram   = "none",
            density.info = "none",
            trace        = "none",
            key          = FALSE,
            margins      = margins,
            lhei         = lhei,
            lwid         = lwid,
            colsep       = 1:p,
            rowsep       = 1:n,
            sepcolor     = sepcolor,
            sepwidth     = sepwidth)
}
