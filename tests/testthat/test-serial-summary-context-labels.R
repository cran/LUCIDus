test_that("serial summary hides Y for non-last stage and labels transition clusters for parallel->early", {
  set.seed(9201)
  n <- 60
  G <- matrix(rnorm(n * 3), nrow = n)
  Z1 <- matrix(rnorm(n * 4), nrow = n)
  Z2 <- matrix(rnorm(n * 4), nrow = n)
  Z3 <- matrix(rnorm(n * 5), nrow = n)
  Y <- rnorm(n)

  Z <- list(list(Z1, Z2), Z3)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "serial",
      G = G,
      Z = Z,
      Y = Y,
      family = "normal",
      K = list(list(2, 2), 2),
      max_itr = 10,
      max_tot.itr = 30,
      tol = 1e-2,
      seed = 9201
    )
  )))

  s <- summary_lucid(fit)
  txt <- capture.output(print(s))

  i_stage1 <- grep("^--- Stage 1", txt)[1]
  i_stage2 <- grep("^--- Stage 2", txt)[1]
  stage1_txt <- txt[i_stage1:(i_stage2 - 1)]
  stage2_txt <- txt[i_stage2:length(txt)]

  expect_false(any(grepl("^\\(1\\) Y", stage1_txt)))
  expect_true(any(grepl("^\\(1\\) Y", stage2_txt)))

  expect_true(any(grepl("previous serial stage", txt, fixed = TRUE)))
  expect_true(any(grepl("Stage1.Layer1.cluster2", txt, fixed = TRUE)))
  expect_true(any(grepl("Stage1.Layer2.cluster2", txt, fixed = TRUE)))
})

test_that("serial summary transition labels are stage-cluster names for early->early", {
  set.seed(9202)
  n <- 60
  G <- matrix(rnorm(n * 3), nrow = n)
  Z1 <- matrix(rnorm(n * 5), nrow = n)
  Z2 <- matrix(rnorm(n * 5), nrow = n)
  Y <- rnorm(n)

  Z <- list(Z1, Z2)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "serial",
      G = G,
      Z = Z,
      Y = Y,
      family = "normal",
      K = list(2, 2),
      max_itr = 10,
      max_tot.itr = 30,
      tol = 1e-2,
      seed = 9202
    )
  )))

  s <- summary_lucid(fit)
  txt <- capture.output(print(s))

  expect_true(any(grepl("Stage1.cluster2", txt, fixed = TRUE)))
})
