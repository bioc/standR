read_decon_fixture <- function(rmNegProbe = FALSE) {
  suppressMessages(readGeoMx(
    test_path("fixtures", "geomx-counts.tsv"),
    test_path("fixtures", "geomx-samples.tsv"),
    rmNegProbe = rmNegProbe
  ))
}

test_that("the default preserves selected log2 expression on the linear scale", {
  spe <- read_decon_fixture()
  # Deliberately differ from raw counts to represent a processed assay.
  linear <- matrix(c(100, 200, 4, 300, 400, 8), nrow = 3,
                   dimnames = dimnames(spe))
  SummarizedExperiment::assay(spe, "logcounts") <- log2(linear)

  out <- prepareSpatialDecon(spe)

  expect_named(out, c("normCount", "backGround"))
  expect_equal(out$normCount, linear)
  expect_equal(unname(out$backGround), matrix(c(4, 4, 4, 8, 8, 8), nrow = 3))
  expect_equal(SummarizedExperiment::assay(spe, "logcounts"), log2(linear))
})

test_that("linear and custom log2 assays are supported explicitly", {
  spe <- geomxNorm(read_decon_fixture(), method = "TMM", log = FALSE)
  linear <- SummarizedExperiment::assay(spe, "logcounts")

  out <- prepareSpatialDecon(spe, isLog = FALSE)
  expect_equal(out$normCount, linear)
  expect_equal(as.numeric(out$backGround[1, ]), as.numeric(linear["NegProbe-WTX", ]))

  SummarizedExperiment::assay(spe, "corrected") <- log2(linear)
  custom <- prepareSpatialDecon(spe, assay2use = "corrected", isLog = TRUE)
  expect_equal(custom, out)
})

test_that("recomputed CPM uses raw counts with matched background scaling", {
  spe <- read_decon_fixture()
  SummarizedExperiment::assay(spe, "logcounts") <-
    SummarizedExperiment::assay(spe, "logcounts") + 10

  out <- prepareSpatialDecon(spe, recomputeCPM = TRUE)
  expected <- matrix(c(10, 20, 1) / 31, nrow = 3) * 1e6
  expected <- cbind(expected, c(30, 40, 2) / 72 * 1e6)
  dimnames(expected) <- dimnames(spe)

  expect_equal(out$normCount, expected)
  expect_equal(as.numeric(out$backGround[1, ]), c(1 / 31, 2 / 72) * 1e6)
  expect_equal(out$backGround[1, ], out$backGround[3, ])
  expect_equal(prepareSpatialDecon(spe, assay2use = "unused", recomputeCPM = TRUE), out)
})

test_that("metadata probes share CPM library sizes and follow ROI order", {
  retained <- prepareSpatialDecon(read_decon_fixture(), recomputeCPM = TRUE)
  spe <- read_decon_fixture(rmNegProbe = TRUE)
  spe <- spe[, c("ROI2", "ROI1")]

  out <- prepareSpatialDecon(spe, recomputeCPM = TRUE)
  genes <- rownames(spe)
  rois <- colnames(spe)
  expect_equal(out$normCount, retained$normCount[genes, rois, drop = FALSE])
  expect_equal(out$backGround, retained$backGround[genes, rois, drop = FALSE])
})

test_that("raw counts can use raw metadata probes without a log conversion", {
  spe <- read_decon_fixture(rmNegProbe = TRUE)
  out <- prepareSpatialDecon(spe, assay2use = "counts")

  expect_equal(out$normCount, as.matrix(SummarizedExperiment::assay(spe, "counts")))
  expect_equal(unname(out$backGround), matrix(c(1, 1, 2, 2), nrow = 2))
})

test_that("raw metadata probes cannot silently supply a processed background", {
  spe <- read_decon_fixture(rmNegProbe = TRUE)
  expect_error(prepareSpatialDecon(spe), "recomputeCPM=TRUE", fixed = TRUE)

  spe <- geomxNorm(spe, method = "TMM", log = FALSE)
  expect_error(prepareSpatialDecon(spe, isLog = FALSE), "recomputeCPM=TRUE", fixed = TRUE)
})

test_that("multiple probe pools retain their own linear background", {
  counts <- matrix(c(10, 1, 20, 2, 30, 3, 40, 4), nrow = 4,
                   dimnames = list(c("GeneA", "NegA", "GeneB", "NegB"), c("ROI1", "ROI2")))
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts, logcounts = log2(counts))
  )

  out <- prepareSpatialDecon(spe, negProbeName = c("NegA", "NegB"),
                            pool = c("A", "A", "B", "B"))
  expect_equal(out$normCount, counts)
  expect_equal(unname(out$backGround), matrix(c(1, 1, 2, 2, 3, 3, 4, 4), nrow = 4))

  expect_error(prepareSpatialDecon(spe, negProbeName = c("NegA", "NegB"), pool = "A"),
               "length(pool)", fixed = TRUE)
})

test_that("missing probes and invalid flags report errors", {
  spe <- read_decon_fixture()
  expect_error(prepareSpatialDecon(spe, negProbeName = "MissingProbe"),
               "must be included", fixed = TRUE)
  expect_error(prepareSpatialDecon(spe, recomputeCPM = NA))
  expect_error(prepareSpatialDecon(spe, recomputeCPM = c(TRUE, FALSE)))
  expect_error(prepareSpatialDecon(spe, isLog = NA))
})
