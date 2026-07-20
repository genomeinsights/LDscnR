# Regression tests for interpolate_cM(), the genetic-map interpolation
# helper used by ld_prune_and_eMLG()'s cM-based run splitting.

build_genetic_map <- function() {
  data.table::data.table(
    Chr = rep(c("Chr1", "Chr2"), each = 4),
    Pos = c(0, 1e6, 2e6, 3e6,   0, 5e5, 1e6, 1.5e6),
    cM  = c(0, 5,   20,  22,    0, 1,   1,   4)
  )
}

test_that("interpolate_cM() linearly interpolates within a chromosome's mapped range", {
  gm <- build_genetic_map()

  ## exact map points should reproduce their own cM value
  out <- interpolate_cM(gm, Chr = rep("Chr1", 4), Pos = c(0, 1e6, 2e6, 3e6))
  expect_equal(out, c(0, 5, 20, 22))

  ## midpoint between two known points -> linear interpolation
  out_mid <- interpolate_cM(gm, Chr = "Chr1", Pos = 5e5)
  expect_equal(out_mid, 2.5)  # halfway between (0, 0) and (1e6, 5)
})

test_that("interpolate_cM() keeps chromosomes independent", {
  gm <- build_genetic_map()

  out <- interpolate_cM(gm, Chr = c("Chr1", "Chr2"), Pos = c(1e6, 1e6))
  expect_equal(out, c(5, 1))  # same Pos, different chromosomes -> different cM
})

test_that("interpolate_cM() clamps out-of-range positions to the nearest endpoint", {
  gm <- build_genetic_map()

  out <- interpolate_cM(gm, Chr = c("Chr1", "Chr1"), Pos = c(-100, 5e6))
  expect_equal(out, c(0, 22))  # clamped to the first/last map point's cM
})

test_that("interpolate_cM() returns NA for chromosomes absent from the map or with <2 points", {
  gm <- build_genetic_map()
  gm_sparse <- rbind(gm, data.table::data.table(Chr = "Chr3", Pos = 100, cM = 0))

  out <- interpolate_cM(gm, Chr = c("Chr1", "Chr9"), Pos = c(1e6, 1e6))
  expect_equal(out[1], 5)
  expect_true(is.na(out[2]))  # Chr9 not in the map at all

  out_sparse <- interpolate_cM(gm_sparse, Chr = "Chr3", Pos = 100)
  expect_true(is.na(out_sparse))  # Chr3 has only 1 point, can't interpolate
})
