# LUCIDus 3.2.1

A correctness follow-up to 3.2.0, focused on `boot_lucid()` and `summary()`.

**Bootstrap CIs now align cluster labels.** Each bootstrap replicate is an
independent refit whose latent cluster *k* need not be the point estimate's
cluster *k*. Replicate coefficient vectors were stacked by name/position anyway,
so any coefficient whose canonical cluster ordering was weakly identified -- a
parallel layer with a near-null cluster-to-outcome effect, an upstream serial
stage ordered by a single omics feature -- could get a bimodal or sign-split
interval spanning orders of magnitude. `boot_lucid()` now permutes every
replicate's clusters to match the reference fit (by posterior-probability
overlap on the resampled subjects) before extracting coefficients, so each
interval brackets its own point estimate. The point estimate (`bootstrap$t0`)
is unchanged.

**`summary()` bootstrap CI tables.**

* The exposure/transition (`(3) E`) and binary-outcome (`(1) Y`) CI tables now
  carry an odds-ratio view -- `OR`, `OR_lower`, `OR_upper` = `exp()` of the
  coefficient-scale estimate and interval bounds -- alongside the
  coefficient-scale columns. `sig` is still computed on the coefficient scale
  (null = 0).
* Headers match the printed scale: the no-CI table says "odds ratio" (an `OR`
  column is shown), the CI table says the columns are coefficients on the
  log-odds scale with a bootstrap interval.
* The `(3) E` CI table now shows the intercept and every covariate (`CoG`) row,
  so it has the same coefficient rows as `summary(fit)` without a bootstrap --
  previously those rows silently vanished when `boot.se` was supplied. The
  covariate row is also now shown in the no-CI `(3) E` table for the parallel
  and serial models, matching the early model.
* `f.binary.early` joins the bootstrap CI to the point estimate by row name
  rather than by position.

**Log-likelihood labeling.** The "Finished LUCID ... model: penalized
log-likelihood = X" completion message printed the *unpenalized* observed-data
value; it now says "observed-data log-likelihood". `em_control$loglik_trace` is
the unpenalized observed-data log-likelihood for every model type (it was the
penalized objective for a penalized early fit), so its last value equals
`$likelihood`; the penalized objective still drives the internal
convergence/dip check.

# LUCIDus 3.2.0

This release follows directly on 3.1.0's correctness pass. It narrows the
exported API to what a user actually calls, adds four extractor functions
for a fitted model's output, marks statistically significant bootstrap CI
rows in `summary()`, fixes several consistency and reliability issues across
model types, and reworks the tutorial vignettes into a polished, user-facing
set built entirely on the public API.

**Breaking change.** Five numerical-stability helpers
(`check_and_stabilize_sigma()`, `safe_solve()`, `safe_log_sum_exp()`,
`check_convergence()`, `safe_normalize()`) and `fill_data()` are no longer
exported. They were internal EM building blocks that happened to have a
CRAN manual page; their behavior is unchanged, and they remain reachable as
`LUCIDus:::fn(...)` for anyone who was calling them directly. The bare
`plot()` S3 generic LUCIDus declared is also removed -- `graphics::plot()`
already provides it, and `plot.early_lucid()` etc. still dispatch normally.

**New: model-type autodetection.** `predict_lucid()` and `boot_lucid()` no
longer require `lucid_model` -- it is now optional and is auto-detected from
`class(model)`. Explicitly passing it still works exactly as before
(including the existing mismatch error when it disagrees with `model`'s
class), so no existing call site is affected; it simply becomes unnecessary
to specify.

**New: extractor functions.** Four new functions read directly off a fitted
model and auto-detect early/parallel/serial themselves, so no branching on
model type is needed at the call site:

* `get_selected_G()` -- the selected exposures, as a named logical vector.
* `get_selected_Z()` -- the selected omics features, shaped by model type
  (a vector, or a list by layer/stage).
* `get_cluster_assignment()` -- each subject's hard cluster assignment.
* `get_top_omics_features()` -- the top-`n` most important omics features
  per layer/stage, ranked by the same criterion `plot_cluster_omic_profile()`
  uses.

**New: significance markers in `summary()`.** Every bootstrap CI table
`summary()` prints when `boot.se`/`se` is supplied -- the exposure-to-cluster
table, the cluster-to-outcome table, and the cluster-specific omics-mean
table -- now includes a `sig` column marking rows whose normal-theory
interval excludes 0 with `"*"`.

**Bug fixes.** The serial model's returned object is now consistent with
early and parallel: it reports a `likelihood` and a `select` field the way
the other two model types always have, a stray unused `Rho_G` penalty is no
longer silently applied where it has no exposures to act on, and its console
progress messages use the same wording and level of detail the other two
model types use. Several input-validation gaps that only showed up for
specific model types or specific entry points (`estimate_lucid()`,
`tune_lucid()`, `predict_lucid()`, `boot_lucid()`) are closed, so a bad input
is now reported with a clear message regardless of which model type or
function it went through, rather than surfacing later as a confusing
downstream error.

**Reliability.** Convergence checking is now handled by a single, shared
routine for every model type, and the optimizer's log-likelihood is checked
at every iteration to confirm it is actually improving, with a warning if it
ever isn't. A handful of tuning and bootstrap code paths that could
previously abort an entire multi-candidate run over one bad candidate now
skip that candidate and continue, matching how the other paths already
behaved.

**Internal.** Every internal (non-user-facing) function in the package now
has real documentation explaining what it does, and source files have been
reorganized and renamed so that related functionality lives together and
file names describe their contents -- neither change affects behavior, and
both are aimed at making the package easier to review and maintain. Several
long-unused, dead code paths were also removed.

**Vignettes.** `helix_early_parallel_workflow.Rmd` is replaced by two
files, `lucid_3models_normal_outcome.Rmd` and
`lucid_3models_binary_outcome.Rmd`, split by outcome family (both still
demonstrate all three model architectures). Every vignette is rewritten as a
user-facing tutorial: custom code that reimplemented feature-selection
extraction, hard cluster assignment, or omics-importance ranking is replaced
by the new extractor functions and `lucid_model` autodetection, sections that
walked through internal-only functions (`fill_data()`, the numerical-stability
helpers) are removed, and each major step gains an explanation of what it
does, why it's needed, and how to read its output.

**Packaging.** The `README.md` is no longer shipped inside the package
tarball (it remains in the source repository), and the paper-reproduction
vignette (`rjournal_paper_examples.Rmd`) is no longer part of the package.
Example, vignette and test runtime is trimmed so a full `R CMD check`
completes well within CRAN's time budget: heavy examples run on small data
subsets with capped EM iterations, vignette bootstrap sizes are reduced, and
essentially the whole model-fitting regression suite is marked
`skip_on_cran()` -- those tests still run in full locally and in CI; only the
fast oracle, input-validation and tightly-capped smoke tests run on CRAN.

# LUCIDus 3.1.0

This release is a correctness pass over the estimation, missing-data,
prediction, model-selection and inference paths, validated against the
statistical definitions in the two LUCID papers:

* Jia, Zhao, Goodrich, Conti (2024). *LUCIDus: An R Package For Implementing
  LUCID With Phenotypic Traits*. The R Journal 16/2. (RJ-2024-012)
* Zhao, Jia, Goodrich, Darst, Conti (2024). *An extension of LUCID incorporating
  incomplete omics data*. Bioinformatics Advances, vbae123.

## Corrections to results reported by earlier versions

These change numbers that version 3.0.x produced. Analyses run with 3.0.x may
need to be re-run.

* **Reported log-likelihood was wrong.** `summary()` reported a value that was
  not the observed-data log-likelihood of the fitted model. For the worked
  example in the R Journal paper (`fit1` on `simulated_HELIX_data`), 3.0.3
  reported `-6721.801`; the correct value at 3.0.3's own maximum-likelihood
  estimates is `-6575.88`. The BIC printed alongside it inherited the error, and
  also counted `K*M^2` covariance parameters rather than the `K*M*(M+1)/2` of
  Eq 13. Both are now correct.

* **Outcome coefficients were mislabelled.** `res_Gamma$beta` held the *absolute*
  cluster means, but `summary()` printed them as if reference-coded
  (`cluster1 = 0`, `cluster2 = <value>`), so the second row read as a
  between-cluster contrast when it was that cluster's own mean. For the paper's
  `fit1`, the printed `0.504` was cluster 2's mean; the actual contrast between
  clusters is `1.044`. Outcome parameters are now stored in both forms and
  printed as an intercept (cluster 1's level) plus explicit contrasts.

* **`predict_lucid()` returned 0-based cluster labels.** `pred.x` was
  `0, ..., K-1` while Eq 21, `summary()`, `inclusion.p` columns and the `mu` /
  `beta` row names all use `1, ..., K`. It now returns `1, ..., K` for early,
  parallel, serial and g-computation. Code that added 1 to `pred.x` must drop
  that adjustment.

* **Penalized BIC now implements Eq 18.** `BICp = -2 logL + (D - DG - DZ) log N`
  subtracts one parameter per *deselected variable*. The previous
  implementation rebuilt the model dimension from only the selected features,
  removing whole mean, covariance and regression blocks, which could select a
  different penalty and therefore a different final model. Note that `D` counts
  the multinomial-logit intercepts, making it larger than Eq 13 as printed by
  `K - 1`; the intercepts are genuinely estimated parameters.

## Bug fixes

* **`K >= 2` was not checked on every path that fits or tunes a model.**
  `lucid()` checked it for early, `tune_lucid()`'s parallel branch checked it,
  and `estimate_lucid()`'s serial branch checked it, but `est_lucid()` itself
  (early and parallel) and `tune_lucid()`'s own early branch did not, so
  calling `estimate_lucid()`/`tune_lucid()` directly with `K = 1` -- both are
  documented as directly callable -- reached deep EM machinery before failing
  with an obscure error instead of a clear one. All four paths now share the
  same check.

* **`predict_lucid()` validated its two branches asymmetrically.** Early
  checked that `model` was an `early_lucid` object; parallel had no analogous
  check. Parallel checked that `Z` had the right number of layers; early
  never checked that `Z` had the right number of columns, so a wrong-width Z
  failed inside `Estep_early()` with a linear-algebra error instead of a
  clear one. Both branches now validate both things.

* **`boot_lucid()` skipped validation the fitting functions it wraps already
  enforce.** None of its three branches called the same G/Y/CoG/CoY
  completeness check `lucid()`/`estimate_lucid()` run, so a missing value
  failed confusingly inside `boot::boot()`'s resampling callback instead of
  with the package's normal top-level error; nor did any branch check that
  `model`'s class matched the `lucid_model` argument. Both checks are now
  applied uniformly across all three branches.

* **`tune_lucid()`'s serial branch could crash on one bad candidate, and had
  an unguarded tie-break.** Early and parallel wrap each grid candidate in
  `try()`, record `NA` on failure, and only `stop()` if every candidate
  failed; serial did neither, so one non-converging `K` (e.g. too many
  clusters for the sample size) crashed the whole tuning run, and picking the
  best model (`which(bic == min(bic))`, no tie-break) would error on a
  genuine BIC tie the way early/parallel are already written to guard
  against. Serial now matches early/parallel's behavior exactly.

* **Serial `Rho_G` was applied to a stage's incoming cluster probabilities as
  if they were exposures.** From stage 2 onward, a serial fit's "G" input is
  the *previous* stage's posterior cluster-membership probabilities, not the
  cohort's actual exposures -- there is nothing there for an exposure-selection
  penalty to select among. `Rho_G` was only reset to 0 when that stage's
  incoming pseudo-G happened to have fewer than 2 columns, so a stage with 2 or
  more forwarded columns (e.g. following a stage with `K >= 4`) was silently
  penalized on a quantity that isn't an exposure. `Rho_G` is now forced to 0 for
  every stage after the first, unconditionally. This changes numeric output for
  that specific combination (nonzero `Rho_G`, a later stage with `K >= 4`
  feeding forward `>= 2` non-reference columns); it does not affect the more
  common case already covered by the old guard.

* **A serial fit's `likelihood` and `select` were commented-out dead code.**
  `R/EM_all.R`'s serial results list referenced `loglik_update`, `selectG`,
  `selectZ` and `Estep_r`, none of which are ever assigned in the serial
  branch -- copy-paste leftovers from the early/parallel branches, left
  commented out. `fit$likelihood` and `fit$select` now exist on a serial fit,
  matching early and parallel: `likelihood` is the sum of each stage's own
  log-likelihood (matching the existing, previously `summary()`-only
  `cal_loglik_serial()`), and `select` is stage 1's own selection (the only
  stage whose `selectG` is a real exposure-selection result, per the fix
  above); the complete per-stage record remains at
  `fit$submodel[[i]]$select` and `fit$submodel[[i]]$likelihood`. The dead `#z
  =` line is removed rather than resurrected: nothing in the package reads
  `z` even where it already exists, on parallel fits.

* **`check_convergence()` was exported, documented and tested, but never
  called.** Early and parallel each ran their own `abs(old - new) < tol`
  convergence check instead. Both now call the shared, tested function.

* **No log-likelihood trace was ever checked for monotonicity.** A decrease
  beyond a small, majorization/selection-aware slack now raises a `warning()`
  naming the iteration (or serial stage) and the two values, for all three
  model types -- previously nothing would have surfaced a real ascent
  violation short of manually inspecting `em_control$loglik_trace`.

* **Serial's verbose output mixed "Sub Model" and "Stage" naming for the same
  concept, and never showed a log-likelihood.** The pre-stage banner now says
  "Stage N/M" consistently, and both the verbose and concise per-stage
  finish lines report that stage's log-likelihood, matching early/parallel's
  own per-iteration reporting.

* **A second round of console-output inconsistencies, found while auditing
  the rest of the package for the same pattern.** `lucid()`'s tuned-selection
  messages and `tune_lucid()`'s `"Tuning LUCID model"` banner printed
  unconditionally instead of respecting `verbose_tune`; `predict_lucid()`'s
  g-computation notice printed unconditionally instead of respecting
  `verbose`, and its `verbose` argument otherwise did nothing for early/
  parallel prediction while serial always printed regardless of it -- early
  and parallel prediction now report progress under the same flag serial
  already did. `predict_lucid()`'s serial per-stage message also still said
  `"Sub Model Number = N"`, the phrasing already replaced with `"Stage N/M"`
  everywhere else. The three "fit completed" messages (`"Success: Early
  Integration LUCID converges!"` / `"Success: LUCID in parallel converges!"`
  / `"Success: LUCID in Serial Model is constructed!"`) and `estimate_lucid()`'s
  own two differently-worded serial "Fitting..." banners are now one
  consistent pattern each.

* **A serial fit's inputs were not validated for completeness at the top
  level.** `check_complete_input()` (the G/Y/CoG/CoY missing-value guard) was
  called from `lucid()` and from the early path, but not from `EM_all.R`'s
  serial orchestration -- a missing value was only caught once stage 1's
  nested fit re-validated it, reporting the error as coming from inside a
  sub-fit. It is now called once at the top of the serial branch.

* **`predict_lucid()` failed whenever `lucid_model` was not named explicitly.**
  The model type was resolved with `match.arg()` only inside the branch tests,
  while the unresolved default -- the full vector
  `c("early", "parallel", "serial")` -- was forwarded to the internal worker,
  whose own choices are `c("early", "parallel")`. Every such call died with
  `'arg' must be of length 1`. This is exactly how the R Journal paper documents
  prediction in Section 3.5 (`predict_lucid(model = fit1, G = ..., Z = ...)`),
  so that example had never run as printed, in this version or in 3.0.3. The
  argument is now resolved once, at entry.

* **Prediction with sporadically missing omics.** Rows containing `NA` were
  passed to a full multivariate Gaussian density. For the early model this
  returned a uniform posterior for every affected row; for the parallel model it
  raised `NA/NaN/Inf in foreign function call`. Both now evaluate the correct
  marginal density on the observed coordinates, which also makes a fully
  unobserved row reduce to the list-wise partition of vbae123 Eq 11.

* **Parallel penalized exposure model discarded posterior mass.** `glmnet`'s
  multinomial fit returns one coefficient block per class, and `f_GtoX()` then
  prepended a further zero reference row, producing `K + 1` states. The E-step
  used only the first `K`, silently dropping roughly 33% of the probability mass
  at `K = 2` and 25% at `K = 3`. Coefficients are now rebased onto a single
  reference class and stored as `K - 1` rows.

* **Parallel penalized omics selection never deselected anything.** The
  indicator OR-ed in `colSums(abs(Sigma)) > 0.001`, which is always true for a
  valid covariance. Selection is now driven by the cluster means, and parallel
  penalized estimation uses the same estimator as the early model (L1-penalized
  means via the graphical-lasso precision matrix).

* **Variable selection could abort before returning a model.** When exactly one
  exposure survived selection, `lucid()` carried the tuned `Rho_G` into the
  refit, where the lasso requires at least two exposures, and the call failed.
  The refit on selected features is now unpenalized; the tuned penalties and the
  original selection are retained on the fitted object (`$Rho`, `$selection`).

* **`init_impute = "lod"` fabricated omics values.** Rows with no measured omics
  at all were filled with `LOD/sqrt(2)` and returned in `fit$Z`, unlike the
  `"mix"` path which restores them. Estimation was unaffected, but the returned
  data were not. All initializers now restore fully-missing rows to `NA`.

* **Bootstrap `t0` was a random refit.** The observed-data statistic re-estimated
  the model under a random seed, so `bootstrap$t0` need not equal the model
  passed in and `summary(fit, boot.se = ...)` could print an `estimate` column
  inconsistent with `summary(fit)`. `t0` is now taken from the supplied model.

* **Failed bootstrap replicates were silent.** Replicates returning non-finite
  estimates were dropped without notice, quietly reducing the effective `R` --
  most likely with missing data, where a resample can contain too few rows with
  complete omics. Failures are now counted and warned about. A small number of
  replicates also raises a warning that the limits are unstable (Davison &
  Hinkley suggest R >= 200 for normal and R >= 800 for percentile intervals).
  Neither suppresses the limits: they are reported as `NA` only when fewer than
  `min_valid` replicates survive, which defaults to 2 -- the point below which
  a standard deviation and an order-statistic interval are undefined rather than
  merely imprecise. `boot_lucid()` gains a `min_valid` argument for callers who
  want a stricter rule.

* **Cross-layer missingness summary.** `check_na(..., "parallel")` reported the
  unweighted mean of per-layer missingness rates rather than the proportion of
  missing cells; a fully-missing 1-feature layer beside a complete 9-feature
  layer reported `0.50` instead of `0.10`.

* **Bootstrap replicate failures no longer print raw error dumps.** The
  per-replicate refit was wrapped in `capture.output()`, which redirects stdout
  while `try()` writes to stderr, so a failed replicate printed a raw
  `Error in est_lucid(...)` trace even though the failure was handled. Two
  sibling call sites already passed `silent = TRUE`; the third now does too.

* **Recovered failures no longer print raw error dumps.** When the penalized
  exposure model could not be fitted (for example `glmnet` reporting
  `Null probability for class 1 < 1.0e-5` on a near-degenerate cluster), the
  code correctly fell back to unpenalized estimation and warned -- but the
  underlying `try()` was not silent, so the raw C++ error text was also printed.
  A single document could emit a dozen of these, each looking like a failure
  when nothing had actually gone wrong. The fallback warning is unchanged; only
  the spurious error text is gone. Affects `Mstep_G()`, the parallel
  `Mstep_GtoX()`, and the early tuning loop.

* Bootstrap confidence-interval tables no longer collapse to a vector when there
  is a single exposure or `K = 2`.

* Penalty tuning no longer performs recursive list indexing when several
  candidate models tie on BIC.

* `glasso` is no longer called when `Rho_Z_Cov = 0`, where it reduces to the
  empirical covariance but emitted a rank warning on every cluster, layer and
  iteration (hundreds of identical warnings per fit).

* **`lucid_serial` fits now report whether they converged.** A serial fit runs
  no EM loop of its own, so its `em_control` carried only the stopping controls
  (`tol`, `max_itr`, `max_tot.itr`). A stage that exhausted `max_itr` was
  therefore invisible unless the user inspected every `submodel[[i]]$em_control`
  by hand. `em_control` now also reports `converged` (`TRUE` only when every
  stage converged), `n_iter` (the total across stages), and the per-stage
  `submodel_converged` and `submodel_n_iter`.

* **`summary()` on a serial fit no longer returns duplicate components.** The
  returned list appended a second copy of `BIC` and `loglik` under names it
  already used, so it held two entries for each. `$BIC` resolved to the first,
  making the copies unreachable, but `str()` output looked malformed and
  name-based iteration visited them twice. The duplicates are removed; `$BIC`,
  `$loglik` and the legacy `$summary.list` alias are unchanged.

* **Missing or non-finite values outside the omics data are now rejected.**
  LUCID models missingness in `Z`; a missing value in `G`, `Y`, `CoG` or `CoY`
  was handled by no code path at all. A single `NA` in `Y` did not error --
  the fit ran on without converging, and where it did return, the outcome
  parameters were badly wrong with nothing in the output to indicate it. These
  four inputs are now checked for `NA`, `NaN` and infinite values before
  fitting, with an error naming the input and the number of offending values.
  `Z` is unaffected and may still be missing.

* **`init_omic.data.model = NULL` now works for the parallel model.** The
  documented automatic-selection value could not be used: `R/em.R` built the
  per-layer geometry names with `rep(init_omic.data.model, length(K))`, and
  `rep(NULL, n)` is `NULL`, so the names were lost between initialization and
  the M-step. `mclust::mstep()` then received `modelName = NULL` and failed,
  reaching the user as "LUCID model fails to converge given current tuning
  parameters". The selected geometry is now read back from the `mclust` fit and
  reported on `init_omic.data.model`, one entry per layer. Random
  initialization, which has no fit to read, falls back to `"EEV"` as the early
  path does. The same line also used `rep()` where `rep_len()` was meant, so a
  per-layer vector of names was expanded to `length(K)^2` entries.

* **`init_impute` and `init_par` are recorded as resolved values.** A serial fit
  taking the defaults stored the unevaluated `c("mix", "lod")` and
  `c("mclust", "random")` -- a record of a choice that was never made, which a
  bootstrap refit reading those fields would carry forward.

* **`predict_lucid()` now enforces one rule for the omics input.** `Z` is
  required for every model type; only `g_computation = TRUE` relaxes it, because
  that mode predicts from the exposure path alone. Previously each model type
  validated separately, and serial got it wrong: the check comparing
  `length(Z)` with `length(K)` ran before the `is.null(Z)` check, and
  `length(NULL)` is 0, so a forgotten `Z` was reported as
  "Z and K should be two lists of the same length" -- sending the user to
  inspect their stage topology rather than the argument they omitted. The
  `is.null` branch was unreachable. All three model types now give the same
  message, and it names `g_computation` so the error states the rule.

* **A single-stage serial model is now declined at prediction with an
  explanation.** `lucid(..., lucid_model = "serial", K = list(2))` fits, but
  `predict_lucid()` failed on it with `object 'pred.y' not found`: the stage
  loop branches on first / middle / last, and a single stage is simultaneously
  first and last, so it takes the first branch and never reaches the assignment
  of `pred.y`. Such a model is a fully equivalent early or parallel model, so
  prediction now stops with a message saying so. *Known limitation*: making
  one-stage serial genuinely predictable requires restructuring the stage loop
  so its inputs depend on being the first stage and its outputs on being the
  last, independently; that is not done here.

## New features

* **`plot_cluster_omic_profile()` shows what the clusters are.** `plot()` draws
  the path structure; `summary()` prints the cluster means as a table. Neither
  answers the question a reader asks first -- which omics features separate the
  clusters, and in which direction. The new function draws that, for early,
  parallel and serial fits, as a cluster-by-feature heatmap in the style used
  for single-cell markers, or as bars. A parallel or serial fit returns one
  plot per omics layer, named, so no figure has to hold every layer at once.

  Features are ranked by an importance measure, `"separation"` by default: the
  spread of the cluster means divided by the typical within-cluster spread. The
  simpler alternatives `"range"` and `"sd"` use the means alone and so cannot
  distinguish a feature that separates the clusters from one that is merely
  noisy. Each returned plot carries the data it drew -- feature, cluster, mean,
  within-cluster SD and score -- as the attribute `"profile_data"`.

  This adds `ggplot2`, `scales` and `grDevices` to Imports.



* **Convergence diagnostics.** Fitted objects carry
  `em_control$converged`, `$n_iter`, `$n_restart` and `$loglik_trace`. Reaching
  `max_itr` without meeting the tolerance previously returned silently as though
  it had converged; it now warns and records `converged = FALSE`.

* **Multiple random starts.** `n_starts` (default `1`) runs independent EM
  initializations and keeps the fit with the highest observed-data
  log-likelihood, with per-start values in `em_control$start_loglik`. EM reaches
  only a local optimum: in testing at `K = 3`, five starts on one dataset spanned
  23 log-likelihood units, and a single start returned a solution whose reported
  exposure-outcome effect was five times too small.

* **Deterministic cluster ordering for unsupervised stages.** Upstream stages of
  a serial model have no outcome to order by, so their cluster numbering was
  seed-dependent and the reported transition coefficients were not reproducible
  across runs on identical data. Such stages are now ordered lexicographically by
  their omics cluster means.

* The imputation step now uses responsibilities recomputed at the updated
  parameters, restoring the ascent property of the E-M-I cycle.

## Documentation

The reference documentation was reviewed against the code function by function.
The substantive corrections:

* **`estimate_lucid()` documented one return object where there are three.** A
  single numbered list described `early`, `parallel` and `serial` alike. It
  listed `em_control` twice with contradictory contents, promised a `likelihood`
  that serial did not yet return at the time (see the Bug fixes entry above --
  it does now) and an `N` that early does not, and omitted `missing_summary`,
  `res_Delta` and `z` altogether. The `@return` is now split into components
  common to every fit and components specific to a model type.

* **`lucid()` duplicated that list** rather than referring to it, and the copy
  had drifted. It now cross-references `estimate_lucid()` and documents the
  `selection` component, which records what the tuned penalties dropped and was
  previously undocumented.

* **`tune_lucid()`** now documents the tuning table's columns, its alignment
  with the returned model list, and that its optimum is the penalized fit --
  `lucid()`, not `tune_lucid()`, performs the unpenalized refit, so estimates
  read directly from `best_model` are shrunk.

* **`fill_data()` described the wrong estimator.** It is not an average of
  per-cluster conditional means; it solves a single precision-weighted system
  over the missing coordinates, which reduces to the Gaussian conditional mean
  only when `K = 1`. The weights are also not the supplied `p` but `p`
  reweighted by each cluster's density at the subject's current values.

* Stub descriptions such as "List containing detailed missing data analysis"
  were replaced with the actual components, conventions and edge-case behaviour
  for `summary_lucid()`, `summary.lucid_parallel()`, `summary.lucid_serial()`,
  the `print` methods, `check_na()`, `analyze_missing_pattern()`,
  `check_imputation_quality()`, `safe_impute()` and the stability utilities.

* `predict_lucid()` now states that `pred.x` runs `1, ..., K`; code written
  against the pre-3.1.0 `0, ..., K - 1` labels must drop its adjustment.

* `plot()` documents its appearance options and their defaults, and states that
  the parallel and serial methods raise an error pending implementation.

* `gen_ci()`, `Istep_Z()` and `summarize_missing_stats()` are internal but were
  generating public help pages; they are now `@noRd`.

## Internal

* `log_dmvnorm_observed()` evaluates Gaussian densities on observed coordinates,
  grouping rows by missingness pattern, and is shared by training and prediction.
* Penalized mean/covariance estimation is factored into a single
  `penalized_cluster_block()` used by both the early and parallel M-steps.
* Test suite rebuilt around independent oracles (`helper-oracle.R`) and a
  canonical data generator (`helper-sim.R`) rather than object-shape assertions.
* `tests/testthat.R` had `test_check()` commented out, so `R CMD check` ran none
  of the suite. It is enabled; the suite runs in roughly three minutes.
* `safe_solve()` and `safe_log_sum_exp()` were each defined three times, and
  `safe_normalize()` and `check_and_stabilize_sigma()` twice, across `utils.R`,
  `stability_utils.R` and `00_stability_preload.R`. The bodies agreed, so
  behaviour was unaffected, but the roxygen sat on a copy that collation then
  overwrote. Consolidated to one definition each in `stability_utils.R`.
* The duplication existed to satisfy a source-time alias in `early_estep.R`,
  which forced a collation order. It is now a call-time wrapper, so no
  load-order coupling remains.
* `safe_solve()` gains `threshold` and `safe_normalize()` gains `min_val`, both
  always supported by the internals but dropped by the exported wrappers.
* `reorder_z()` / `reorder_lucid()` in `R/summary.R` were dead code -- never
  called anywhere, and referencing field names (`res_Delta$Delta`,
  `res_Mu_Sigma`) that no current class has. Removed.
* Early's per-iteration verbose message now says `"log-likelihood"` /
  `"penalized log-likelihood"` instead of `"loglike"` / `"penalized loglike"`,
  matching parallel's own wording; this is a console-output-only change, not a
  numeric one.
* A second orphaned cluster-reordering subsystem removed: `get_ref_cluster()`,
  `reorder_Delta()`, `reorder_Mu_Sigma()`, `reorder_Beta()` (`R/summary.R`,
  the functions the already-removed `reorder_z()`/`reorder_lucid()` used to
  call), `print.sumlucid_auxi()` and `print.auxi.serial.scen1/2/3()` (S3
  methods registered in `NAMESPACE` for classes never assigned to any object,
  so they could never dispatch), `parallel_marginal_design()`,
  `make_positive_definite()`/`is_positive_definite()` (superseded by
  `check_and_stabilize_sigma()`/`safe_solve()`), `fill_data_help1()`/
  `fill_data_help2()`, `initialize_Beta()`, and `early_estep()` (a full
  duplicate of `Estep_early()`) -- all confirmed to have zero call sites
  anywhere in the package. `check_K()`, found dead in the same sweep, was kept
  and wired in instead (see Bug fixes); `plot.lucid_parallel()`'s ~100 lines
  of Sankey-plotting code after its unconditional `stop()` (unreachable) were
  also removed, and `plot.lucid_serial`/`plot.lucid_parallel` now share the
  same trimmed shape.
