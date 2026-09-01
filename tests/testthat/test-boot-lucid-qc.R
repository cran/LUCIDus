skip_on_cran()

# boot_lucid() used to skip the validation the fitting functions it wraps
# already enforce: none of its three branches called check_complete_input()
# on G/Y/CoG/CoY, so a missing value failed confusingly inside boot::boot()'s
# resampling callback instead of with the package's normal top-level error.
# It also never checked that `model`'s class matched the `lucid_model`
# argument, unlike pred_lucid()'s (partial) class check.

fit_for_boot <- function(model, n_layer = 2L, seed = 61) {
  d <- sim_lucid(model, N = 200, K = 2, M = 4, n_layer = n_layer,
                sep = 2.5, family = "normal", seed = seed)
  K <- if (model == "early") 2 else as.list(rep(2, n_layer))
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = model,
                family = "normal", K = K, init_omic.data.model = "VVV")
  ))))
  list(d = d, fit = fit)
}

test_that("boot_lucid() rejects a missing value in Y the same way fitting does, for all three model types", {
  for (model in c("early", "parallel", "serial")) {
    fx <- fit_for_boot(model)
    Y_bad <- fx$d$Y
    Y_bad[3] <- NA
    expect_error(
      suppressWarnings(suppressMessages(invisible(capture.output(
        boot_lucid(G = fx$d$G, Z = fx$d$Z, Y = Y_bad, lucid_model = model,
                  model = fx$fit, R = 5)
      )))),
      "Y contains 1 missing value",
      info = model
    )
  }
})

test_that("boot_lucid() rejects a model/lucid_model class mismatch, for all three model types", {
  fx_early <- fit_for_boot("early", n_layer = 1L)
  fx_parallel <- fit_for_boot("parallel")

  # an early_lucid model passed with lucid_model = "parallel"
  expect_error(
    boot_lucid(G = fx_early$d$G, Z = fx_early$d$Z, Y = fx_early$d$Y,
              lucid_model = "parallel", model = fx_early$fit, R = 5),
    "lucid_parallel"
  )
  # a lucid_parallel model passed with lucid_model = "early"
  expect_error(
    boot_lucid(G = fx_parallel$d$G, Z = fx_parallel$d$Z, Y = fx_parallel$d$Y,
              lucid_model = "early", model = fx_parallel$fit, R = 5),
    "early_lucid"
  )
  # a lucid_parallel model passed with lucid_model = "serial"
  expect_error(
    boot_lucid(G = fx_parallel$d$G, Z = fx_parallel$d$Z, Y = fx_parallel$d$Y,
              lucid_model = "serial", model = fx_parallel$fit, R = 5),
    "lucid_serial"
  )
})

test_that("lucid_model is auto-detected from model's class when omitted", {
  for (model in c("early", "parallel", "serial")) {
    fx <- fit_for_boot(model)
    set.seed(1)
    with_arg <- suppressWarnings(boot_lucid(G = fx$d$G, Z = fx$d$Z, Y = fx$d$Y,
                                            lucid_model = model, model = fx$fit, R = 5))
    set.seed(1)
    without_arg <- suppressWarnings(boot_lucid(G = fx$d$G, Z = fx$d$Z, Y = fx$d$Y,
                                               model = fx$fit, R = 5))
    expect_equal(with_arg, without_arg, info = model)
  }
})
