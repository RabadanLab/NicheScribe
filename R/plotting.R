# Plot colors

.colsR21 <- c(
  "CM_MSCL" = "#0CB702",
  "Endo_MSCL" = "palegreen",
  "Endo_MSCS" = "purple",
  "Endo_mMPr" = "orange",
  "Endo_CD24_neg_MPr" = "red",
  "Endo_CD24_pos_MPr" = "darkturquoise",
  "CM_SEC" = "#F8766D",
  "Endo_AEC" = "#ABA300",
  "Endo_SEC" = "navy",
  "CM_AEC" = "brown"
)

.colsR22 <- c(
  "MSCL" = "#0CB702",
  "Adipo_MSCL" = "palegreen",
  "Osteo_MSCL" = "darkgreen",
  "MSCS" = "purple",
  "OPr1" = "#FF61CC",
  "OPr2" = "violetred",
  "OPr3" = "red",
  "CPr1" = "dodgerblue",
  "CPr2" = "mediumturquoise",
  "FPr" = "plum1",
  "Pericyte" = "slateblue",
  "aEC" = "#ABA300",
  "lsEC1" = "#F8766D",
  "hEC2" = "navy",
  "Doublet" = "black",
  "hEC1" = "skyblue",
  "lsEC2" = "pink3",
  "rEC" = "orange"
)

#' NicheScribe UMAP plotting function
#'
#' Plots cells on NicheScribe-generated UMAP with NicheScribe annotations and standard colors
#'
#' @param query Seurat object to plot from - must have NicheScribe annotations and UMAP
#' @param layer Annotation layer to plot. Either `"R2.1"` (uses `hash_label_pred`) or `"R2.2"` (uses `manual_annot_pred`).
#' @param ... Additional arguments passed to Seurat::DimPlot().
#'
#' @return Seurat DimPlot using NicheScribe UMAP reduction, selected annotations, and standard color palette
#'
#' @export
NichePlot <- function(query, layer = c("R2.1", "R2.2"), ...) {

  layer <- match.arg(layer)

  if (layer == "R2.1") {
    group.by <- "hash_label_pred"
    cols <- .colsR21
  } else {
    group.by <- "manual_annot_pred"
    cols <- .colsR22
  }

  if (!"nichescribe_umap" %in% names(query@reductions)) {
    rlang::abort("nichescribe_umap reduction not found.")
  }

  if (!group.by %in% colnames(query@meta.data)) {
    rlang::abort(group.by, " metadata column not found.")
  }

  Seurat::DimPlot(
    query,
    reduction = "nichescribe_umap",
    group.by = group.by,
    cols = cols,
    ...
  ) +
    ggplot2::coord_cartesian(
      xlim = c(-14, 12),
      ylim = c(-15, 15)
    ) +
    ggplot2::scale_x_continuous(breaks = c(-10, -5, 0, 5, 10)) +
    ggplot2::scale_y_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15)) +
    ggplot2::ggtitle(paste0("NicheScribe - ", layer))
}

#' NicheScribe iMSC-L score plotting function
#'
#' Plots the inflammatory MSC-L score on the NicheScribe UMAP.
#' Optionally splits the plot by a metadata variable.
#'
#' @param query Seurat object containing a nichescribe_umap reduction and
#'   iMSCL_score metadata column.
#' @param group Optional metadata column to split the plot by.
#' @param ... Additional arguments passed to Seurat::FeaturePlot().
#'
#' @return A ggplot or patchwork object.
#'
#' @export
InflammationPlot <- function(query, group = NULL, ...) {

  # Common colour scale
  lims <- range(query$iMSCL_score, na.rm = TRUE)

  # No grouping
  if (is.null(group)) {

    return(
      suppressMessages(
        Seurat::FeaturePlot(
          query,
          reduction = "nichescribe_umap",
          features = "iMSCL_score",
          combine = TRUE,
          ...
        ) +
          ggplot2::scale_color_gradientn(
            colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
            limits = lims
          ) +
          ggplot2::coord_cartesian(
            xlim = c(-14, 12),
            ylim = c(-15, 15)
          ) +
          ggplot2::scale_x_continuous(breaks = c(-10, -5, 0, 5, 10)) +
          ggplot2::scale_y_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15)) +
          ggplot2::ggtitle("iMSC-L score")
      )
    )
  }

  # Check metadata column
  if (!group %in% colnames(query@meta.data)) {
    rlang::abort(group, " is not a metadata column.")
  }

  groups <- sort(unique(query[[group]][,1]))

  plots <- lapply(groups, function(g) {

    q <- subset(
      query,
      cells = SeuratObject::Cells(query)[query[[group]][,1] == g]
    )
    suppressMessages(
      Seurat::FeaturePlot(
        q,
        reduction = "nichescribe_umap",
        features = "iMSCL_score",
        combine = FALSE,
        ...
      )[[1]] +
        ggplot2::ggtitle(as.character(g)) +
        ggplot2::scale_color_gradientn(
          colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
          limits = lims
        ) +
        ggplot2::coord_cartesian(
          xlim = c(-14, 12),
          ylim = c(-15, 15)
        ) +
        ggplot2::scale_x_continuous(breaks = c(-10, -5, 0, 5, 10)) +
        ggplot2::scale_y_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15))
    )
  })

  patchwork::wrap_plots(plots) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}
