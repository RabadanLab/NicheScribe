NicheScribeData <- new.env(parent=emptyenv())

.reference_repo <- "RabadanLab/NicheScribe"
.reference_version <- "v1.0.0"

.reference_path <- function() {
  file.path(tools::R_user_dir("NicheScribe", which = "data"), "reference")
}

.onAttach <- function(libname, pkgname) {
  if (!validate_reference_files()) {
    packageStartupMessage(
      "NicheScribe reference data not found or out of date.\n",
      "Run NicheScribe::download_reference() to install required data."
    )
  }
}

read_reference_manifest <- function(path=.reference_path()) {
  mf <- file.path(path, "manifest.json")
  if (!file.exists(mf)) return(NULL)
  jsonlite::read_json(mf, simplifyVector = TRUE)
}

validate_reference_manifest <- function(manifest) {
  required <- c("schema_version", "reference_version", "files")
  if (!all(required %in% names(manifest))) {
    rlang::abort("Invalid NicheScribe reference manifest.")
  }
  if (!identical(manifest$reference_version, .reference_version)) {
    return(FALSE)
  }
  TRUE
}

validate_reference_files <- function(path=.reference_path()) {
  manifest <- read_reference_manifest(path)
  if (is.null(manifest) || !validate_reference_manifest(manifest)) {
    return(FALSE)
  }

  for (i in seq_len(nrow(manifest$files))) {
    filename <- manifest$files$filename[i]
    expected_hash <- manifest$files$sha256[i]
    fp <- file.path(path, filename)

    if (!file.exists(fp)) return(FALSE)

    observed_hash <- digest::digest(fp, algo="sha256", file=TRUE)
    if (!identical(observed_hash, expected_hash)) return(FALSE)
  }
  TRUE
}

#' Download NicheScribe reference data
#'
#' @description Downloads and installs the required reference datasets for
#' NicheScribe from GitHub. File integrity is verified via SHA256 checksums.
#'
#' @param force Logical. If `TRUE`, forces a re-download even if the data is already up-to-date and validated. Defaults to `FALSE`.
#'
#' @return Invisible `TRUE` if successful.
#' @export
download_reference <- function(force=FALSE) {
  destination <- .reference_path()

  if (!force && validate_reference_files(destination)) {
    rlang::inform("NicheScribe reference is already up to date.")
    return(invisible(TRUE))
  }

  if (!force && interactive()) {
    ans <- utils::askYesNo(paste0("Download reference data (~160 MB) to ", destination, "?"))
    if (!isTRUE(ans)) return(invisible())
  }

  tmp <- tempfile("NicheScribe-reference-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive=TRUE), add=TRUE)

  tryCatch({
    # Download and check manifest
    manifest_file <- "manifest.json"
    piggyback::pb_download(
      file=manifest_file,
      repo=.reference_repo,
      tag=.reference_version,
      dest=tmp
    )

    manifest <- jsonlite::read_json(
      file.path(tmp, manifest_file),
      simplifyVector=TRUE
    )

    if (!validate_reference_manifest(manifest)) {
      rlang::abort("Downloaded manifest does not match this package version.")
    }

    # Download files specified in manifest
    for (i in seq_len(nrow(manifest$files))) {
      piggyback::pb_download(
        file=manifest$files$filename[i],
        repo=.reference_repo,
        tag=.reference_version,
        dest=tmp
      )
    }

    # SHA256 verification
    if (!validate_reference_files(tmp)) {
      rlang::abort("Downloaded reference files failed integrity checks.")
    }

    # Copyfiles to destination, overwriting existing files
    dir.create(destination, recursive=TRUE, showWarnings=FALSE)
    files_to_copy <- list.files(tmp, full.names=TRUE)
    file.copy(
      from=files_to_copy,
      to=destination,
      overwrite=TRUE
    )

    rlang::inform(paste("Sucessfully installed NicheScribe reference", .reference_version))
    invisible(TRUE)
  }, error = function(e) {
    rlang::abort(paste("Reference download failed:", conditionMessage(e)))
  })
}

#' Load NicheScribe reference data
#'
#' @description Loads the validated NicheScribe reference datasets into the
#' dedicated package environment (`NicheScribeData`).
#'
#' @return An invisible environment containing the required reference objects.
#' @export
load_data <- function() {
  path <- .reference_path()

  if (!validate_reference_files(path)) {
    rlang::abort("NicheScribe data missing or corrupted. Run NicheScribe::download_reference().")
  }

  # Only load if not already in the environment
  req_vars <- c("R1", "R2", "umap_model", "imscl_genes")
  if (!all(vapply(req_vars, exists, logical(1), envir=NicheScribeData, inherits=FALSE))) {
    manifest <- read_reference_manifest(path)

    for (i in seq_len(nrow(manifest$files))) {
      f <- manifest$files$filename[i]
      fp <- file.path(path, f)
      if (f == "ref.rds") {
        ref_obj <- readRDS(fp)
        assign("ref", ref_obj, envir=NicheScribeData)
      } else if (f == "umap_model.uwot") {
        uwot_obj <- uwot::load_uwot(fp)
        assign("umap_model", uwot_obj, envir=NicheScribeData)
      }
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

    # Cleanup temporary loaded raw object
    rm(list="ref", envir=NicheScribeData)
  }

  invisible(NicheScribeData)
}
