cellwiseheatmap <- function (cellwiseoutliers,
                             data,
                             col = c("blue", "lightgray", "red"),
                             col.scale.factor = 1,
                             notecol.outlier = "white",
                             notecol.clean = "black",
                             notecex = 1,
                             margins = c(9.5, 14),
                             lhei = c(0.5, 15),
                             lwid = c(0.1, 3.5),
                             sepcolor = "white",
                             sepwidth = c(0.01, 0.01)) {

  n <- nrow(cellwiseoutliers)
  p <- ncol(cellwiseoutliers)

  rownames(cellwiseoutliers) <- rownames(data)
  colnames(cellwiseoutliers) <- colnames(data)

  notecol <- ifelse(c(t(cellwiseoutliers[n:1, ])) != 0, notecol.outlier, notecol.clean)
  ncellcols <- length((cellwiseoutliers))

  realpow <- function (x, rad) {
    sign(x) * (abs(x))^rad
  }

  heatmap.2(realpow(as.matrix(cellwiseoutliers), col.scale.factor),
            scale        = "none",
            col          = colorpanel(n = ncellcols, low = col[1], mid = col[2], high = col[3]),
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
