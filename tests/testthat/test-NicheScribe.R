test_that("NicheScribe snapshot test on test data", {
  path <- test_path("testdata", "test_data.rds")
  seu <- readRDS(path)

  seu <- NicheScribe(seu)

  expect_true(all(c(
    "nonhematopoietic_pred",
    "hash_label_pred",
    "manual_annot_pred"
  ) %in% colnames(seu@meta.data)))

  expect_snapshot(list(
    nonhematopoietic = seu$nonhematopoietic_pred,
    hash = seu$hash_label_pred,
    manual = seu$manual_annot_pred
  ))
})
