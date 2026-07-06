NicheScribeData <- new.env(parent=emptyenv())

ref.files <- list(
  ref = "20260702_seurat_ref_slim.rds",
  umap_model = "20260524_umap_model.uwot"
)

.onAttach <- function(libname, pkgname) {
  if (!.data_available()) {
    packageStartupMessage(
      "NicheScribe reference data not found.\n",
      "Run NicheScribe::download_reference() to install required data."
    )
  }
}

.data_available <- function(dataDir = NA) {
  if (is.na(dataDir)) {
    dataDir <- tools::R_user_dir("NicheScribe", which = "data")
  }

  all(file.exists(file.path(dataDir, ref.files)))
}

get_data <- function(dataDir = NA) {
  if (is.na(dataDir)) {
    dataDir <- tools::R_user_dir("NicheScribe", which = "data")
  }

  if (!.data_available(dataDir)) {
    rlang::abort("NicheScribe data not found. Run NicheScribe::download_reference().")
  }

  if (!all(vapply(c("R1", "R2", "umap_model", "imscl_genes"), exists, logical(1), envir=NicheScribeData, inherits=FALSE))) {
    for (fn in names(ref.files)) {
      f <- ref.files[[fn]]
      fp <- file.path(dataDir, f)
      obj <- switch(
        tools::file_ext(f),
        rds = readRDS(fp),
        uwot = uwot::load_uwot(fp),
        rlang::abort(paste("Unknown file extension:", f))
      )
      assign(fn, obj, envir=NicheScribeData)
    }

    NicheScribeData$R1 <- NicheScribeData$ref$ref
    stopifnot(all(SeuratObject::Cells(NicheScribeData$R1) == SeuratObject::Cells(NicheScribeData$ref$bm_pca_model)))
    NicheScribeData$R1[["pca"]] <- NicheScribeData$ref$bm_pca_model
    NicheScribeData$R2 <- subset(NicheScribeData$R1, cells=SeuratObject::Cells(NicheScribeData$R1)[NicheScribeData$R1$nonhematopoietic == "nonhematopoietic"])
    stopifnot(all(SeuratObject::Cells(NicheScribeData$R2) == SeuratObject::Cells(NicheScribeData$ref$stroma_pca_model)))
    suppressWarnings({
      NicheScribeData$R2[["pca"]] <- NicheScribeData$ref$stroma_pca_model
    })
    stopifnot(all(SeuratObject::Cells(NicheScribeData$R2) == rownames(NicheScribeData$umap_model$embedding)))
    NicheScribeData$imscl_genes <- NicheScribeData$ref$imscl_genes
    rm(list="ref", envir=NicheScribeData)
  }

  invisible(NicheScribeData)
}

#' Download reference data
#'
#' Downloads the reference data from the GitHub repo.
#'
#' @param dataDir Specifies where the reference data should be downloaded. Should essentially never be altered by the user directly since the package looks for the data in a specific location.
#' @export
download_reference <- function(dataDir = NA) {
  if (is.na(dataDir)) {
    dataDir <- tools::R_user_dir("NicheScribe", which = "data")
  }

  if (interactive()) {
    ans <- utils::askYesNo(paste0("Download reference data (~160 MB) to ", dataDir, "?"))
    if (!ans) return(invisible())
  }

  if (!dir.exists(dataDir)) {
    dir.create(dataDir, recursive = TRUE)
  }

  tryCatch({
    for (f in ref.files) {
      piggyback::pb_download(
        file=f,
        repo="RabadanLab/NicheScribe",
        dest=dataDir,
        overwrite=TRUE
      )
    }
    rlang::inform("Reference data successfully downloaded.")
  }, error = function(e) {
    rlang::abort(paste("Download failed:", conditionMessage(e)))
  })
}
