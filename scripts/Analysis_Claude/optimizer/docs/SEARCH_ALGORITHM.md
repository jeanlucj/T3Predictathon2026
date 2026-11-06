# How the search works

The complete rules the optimizer follows when it decides what to evaluate next, what counts as
a good configuration, and how much confidence its ranking carries. `DESIGN.md` §3 is the
summary; this is the full statement, written so that a reader who knows no part of it can
follow every step, and so that no sentence has to be checked against the code to be trusted.

**Scope.** This document covers the *search*: the loop over configurations. It does not cover
what a configuration *does* (the six subtasks -- `DESIGN.md` §1), nor the T3 data layer
(`LESSONS.md` §1), nor how a run is launched (`../RUNBOOK_*.md`).

**How it was verified.** Every rule below is quoted from, or derived from, the code at the file
and line given. Every number labelled *measured* comes from the seven optimization runs in
`state/evals_backup.sqlite` as of 2026-10-30 (build 0.8.8), re-deriving the optimizer's own
estimates with the optimizer's own functions over prefixes of each run's evaluation stream.
Where a rule's consequence has not been measured, it says so.

---

## 1. Glossary

Read this once; the rest of the document uses these words in exactly these senses.

| term | meaning |
|---|---|
| **subtask** | One of the six stages of a genomic-prediction pipeline (`SUBTASKS`, `R/config_space.R:26`): `train_select`, `pheno_prep`, `geno_select`, `kernel`, `model`, `predict_post`. |
| **configuration** (config) | One complete choice of method + parameters for all six subtasks, as a flat named list. Parameters that do not apply to the chosen method are `NA`. |
| **`config_hash`** | The identity of a configuration: a hash of its canonical full-precision JSON (`R/config_space.R:196`). Two configurations are "the same" exactly when this matches. |
| **block** | A subtask's method key plus all its parameter keys, moved as a unit by crossover and mutation. |
| **focal trial** | The trial being predicted in one evaluation. |
| **training trials** | The other trials a configuration pulls phenotypes from. Never restricted by the target domain. |
| **scheme** | `CV0` or `CV00`: the cross-validation rule that decides which accessions may inform the prediction. One run optimizes exactly one scheme. |
| **evaluation** | One run of one configuration against one focal trial under one scheme. The unit of cost and the unit of evidence. One row of the `evals` table. |
| **`score`** | An evaluation's result: the Pearson correlation between predicted and observed per-accession BLUEs on the focal trial (`score_predictions`, `R/evaluate.R:110`). `NA` when the evaluation failed. |
| **`n_test`** | The number of accessions that correlation was computed over -- joined, and finite in both prediction and observation (`R/evaluate.R:133`). The precision of a score. |
| **usable** | A stored row that has a finite score **and** `n_test >= min_test_acc` (**10**). Only usable rows enter any estimate, and the same floor decides whether a score is computed at all -- §3 (`.fisher_weights`, `R/optimizer.R:127`). |
| **target domain** | The declared restriction on which trials may be *focal*: any of programs, years, locations, trial names (`settings$target_domain`, applied by `.apply_target_domain`, `R/data_access.R:391`). |
| **universe** | The target domain resolved **once per run** into a fixed set of trial ids and stored (`pin_trial_universe`, `R/data_access.R:441`). Membership is a set test on an immutable id, so later changes to a trial's metadata cannot move a stored row in or out. |
| **run_id** | The hash that identifies a run: the run-defining settings plus the build (`run_id_for`, `R/store.R`). Same settings, same build => same run, resumed, inheriting its pinned trial universe. **It is provenance, not scope**: nothing in selection reads it, so a new `run_id` does not hide the old rows. §13. |
| **target** | Informally, what a run is optimizing: the universe plus the scheme. Both are recoverable from the `runs` table (`universe_json`, `settings_json`). |
| **store** | `state/evals.sqlite`. A **global archive** of every evaluation ever run, across all targets. Each run reads only its own slice. |
| **Fisher z** | `z = atanh(score)`, the variance-stabilizing transform for a correlation. All model fitting happens on this scale; results are mapped back with `tanh`. |
| **weight** (`w`) | `1/w = 1/max(n_test - 3, 1) + 1/weight_cap` (`.eval_weight`, `R/optimizer.R:122`). Two precisions in series: the sampling variance of a Fisher-z-transformed correlation, `1/(n-3)`, plus a floor `1/weight_cap` for config × trial interaction, which does **not** shrink as the trial grows. `weight_cap` is therefore the weight an unboundedly large trial gets, and reads as a trial size — 50 means "no trial counts for more than one with 53 accessions". Governs the aggregate **and** the surrogate. §4.2, `LESSONS.md` #36. |
| **`sd_trial`** | The standard deviation of the per-trial random effect: how much trials differ from each other in difficulty. |
| **`sd_config`** | The standard deviation of the per-configuration random effect: how much the configurations *in the fit* differ from each other in quality. **A property of the population of configurations evaluated, not a constant of nature.** |
| **`sd_resid`** | The residual standard deviation **at unit weight**. To get the noise of one evaluation at a typical trial, divide by `sqrt(median(w))`. Comparing it to the other two as printed is a mistake (`LESSONS.md` #19). |
| **pooled estimate** | A configuration's weighted mean score on the z scale, back-transformed. Used when the random-effects fit cannot run. |
| **BLUP** | Best Linear Unbiased Predictor: the configuration's estimated quality from the random-effects fit, after removing the trial effect and **shrinking** toward the overall mean in proportion to how little evidence there is. This is `mean_score`, the number everything downstream ranks on. |
| **shrinkage** | Pulling an estimate toward the mean because the data behind it are weak. The weight the BLUP puts on a configuration's own data is `lambda` in §4.5. |
| **`se`** | The **posterior standard deviation** of a configuration's random effect, on the z scale (`condVar`, `R/optimizer.R:264`). Not the standard error of a mean. Supplied *only* by the random-effects fit. |
| **UCB** | Upper confidence bound, `tanh(z_mean + contender_z * se)`. An optimistic view of a configuration, formed on the z scale the fit produced and transformed back once. §7. |
| **contender** | A configuration whose UCB still reaches the leader's `mean_score`, i.e. one whose uncertainty still covers the leader's estimate -- capped at 8 (`.contenders`, `R/optimizer.R:402`). |
| **leader** | The configuration with the highest `mean_score` at a given moment. Not a stored object; recomputed from scratch at every step. |
| **incumbent** | The leader restricted to configurations with at least `incumbent_min_reps` usable evaluations. **Report-only** (`R/report.R:119`); the search never consults it. |
| **run** | A history of evaluations that belong together. What "together" means is a setting: under `run_identity = "target"` (the default) a run is one optimization target and its policy may change within it; under `"experiment"` a run is one target *and* one fixed policy. §13. |
| **reference population** | The configurations drawn WITHOUT reference to the surrogate -- the Phase 3 random draws, the Phase 3b persistent draws and the Phase 5 fallback, read off `evals.source`; seeds excluded. A sample of the configuration *space*, as against the full evaluated set, which is a sample of the search's *trajectory*. `sd_config` is estimated on it (§11). |
| **parent** | A configuration used for crossover and mutation. The in-play set of §7, ordered by UCB, capped at `n_parents` (32) and floored at `n_elites` (8) -- `get_parents`, `R/optimizer.R`. No replication requirement. |
| **elite** | One of the top `n_elites` (8) configurations by `mean_score` (`get_elites`). Since 2026-11-04 this is only the **floor** under the parent pool, for the case where the random-effects fit did not run and there is no `se` to rank on; before that it *was* the parent pool (`LESSONS.md` #42). |
| **surrogate** | The model that predicts a configuration's score from its features, so candidates can be ranked without running them. A bag of 300 `rpart` regression trees (`R/surrogate.R`). |
| **EI** | Expected Improvement: the acquisition function that converts the surrogate's (mean, uncertainty) into one number per candidate (`expected_improvement`, `R/surrogate.R:155`). |
| **backlog** | Configurations owed a further evaluation, in two tiers: `base` (below the floor) and `extra` (contenders above the floor) -- `.replication_backlog`, `R/optimizer.R:384`. |
| **churn** | Defined in §9. |

---

## 2. What is being optimized

The objective is **the expected score of a configuration on a trial drawn at random from the
target domain**, under one scheme.

Four consequences follow, and most of the search's behaviour is downstream of them:

1. **The objective is an average over a population of trials**, so a single evaluation is a
   sample, not a measurement. Two configurations evaluated on different trials are not directly
   comparable until the trial effect is removed.
2. **It is expensive.** One evaluation downloads phenotypes and genotypes (cached per trial),
   builds a kernel and fits mixed models.
3. **The configuration space is infinite.** Five of the six subtasks carry continuous
   parameters sampled from a uniform or log-uniform range (`R/config_space.R:26-108`), so the
   search can never be exhaustive and "have we tried everything" is never a question.
4. **The answer is a distribution, not an identification.** Every quantity the search produces
   is an estimate carrying uncertainty, and on this data that uncertainty is large relative to
   the differences being estimated (§4.6). "The best configuration" is never established; what
   is estimated is each configuration's effect and how precisely it is known.

**Read every ranking in this document in that light.** The search selects, reports and exploits
a configuration because its estimated effect is highest *at that moment and under the current
fit* -- a best guess, revised by the next evaluation, never a decision that closes a question.
*Leader*, *elite*, *parent*, *contender* and *incumbent* name positions in a ranking; none of
them asserts that a configuration has been shown to be better than another, and *contender* and
*parent* are explicitly the positions that say the opposite -- a configuration is in play
precisely because it has **not** been ruled out. Three
places carry the uncertainty explicitly: the shrinkage of §4.5 (how much of a configuration's
own data its estimate is willing to use, and how precise that estimate can ever get), the
contender rule of §7 (which configurations' uncertainty still covers the leader's estimate), and
the resolution series of §9 (whether that set is narrowing at all). Everywhere else the search
ranks on a point estimate and the uncertainty is implicit.

Only the *focal* trial is restricted by the target domain; a configuration may pull training
trials from anywhere in the database.

---

## 3. The unit of evidence: one evaluation

`optimizer_step` (`run_optimizer.R:23`) does exactly this, once per iteration:

1. `choose_config` picks a configuration (§5).
2. `choose_trial` picks a focal trial for it (§6).
3. `claim_eval` atomically claims the `(config_hash, trial_id, scheme)` triple. If another
   worker already holds it, pick again -- up to 3 attempts, then give the step up. This is what
   makes concurrency safe; `choose_trial`'s own avoidance of claimed trials is advisory only
   (`LESSONS.md` #31).
4. `evaluate_config_on_trial` runs the pipeline and scores it.
5. `store_eval` writes one row: the config JSON, trial id, scheme, score, `n_test`, status,
   failure reason, timing, memory, build, `run_id`, and **`source`** -- which phase of
   `choose_config` proposed the configuration. `source` is what distinguishes the rows chosen
   *independently of the surrogate* (`seed:*`, `random_init`, `random_fallback`) from the rows
   the surrogate chose, and those are the only unbiased sample of the configuration universe the
   store contains -- see §11 on why that matters for `sd_config`. It has been logged since the
   first build and **stored only since 2026-11-03**, so archived rows carry `NA` and cannot be
   classified (`LESSONS.md` #38).

A row is written **whether or not the evaluation succeeded**. That matters in §5: the
replication floor counts rows, not successes.

For a score to exist at all, at least `min_test_acc` (**10**) accessions must have both a
finite prediction and a finite observation; below that the row is stored with
`status = "too_few_overlap"` and `score = NA` -- `n_test` still records the overlap that was
found, so the shortfall is visible in the store. The same number decides whether a score that
does exist may be used: `.fisher_weights` marks a row `.usable` only at or above it
(`R/optimizer.R:127`), and since 2026-11-03 the surrogate's training rows pass through the same
helper and inherit the same floor (§8.1). One setting, applied everywhere, so the ends cannot
drift apart.

> **Why it is one setting.** It was two until 2026-11-02 -- feasibility at 5 in
> `score_predictions`, usability at 10 hardcoded in `aggregate_scores` -- and the band between
> them was dead weight. **Measured** on run `9a2664962c01`: six of its fourteen trials returned
> `n_test = 7` on every single evaluation, so 2029 of 5631 successful evaluations (36%) cost
> full compute, were stored as `ok` with a real score, and were then invisible to every
> estimate. The maximum attainable replication for any configuration in that run was 8 usable
> evaluations, not 15, and nothing reported it. See `docs/LESSONS.md` #34.
>
> Rows below the floor still exist in the archive, written before the change. `.usable` is what
> keeps them out of the estimates; it is not dead code.

---

## 4. From scores to a ranking

This is `aggregate_scores` (`R/optimizer.R:149`), and everything else in the search consumes
its output.

### 4.1 Which rows are in

The run's slice of the global store, in this order (`choose_config`, `R/optimizer.R:518`):

- `filter_evals_to_universe` -- only trials in this run's pinned universe;
- `filter_evals_to_scheme` -- only this run's `optimize_scheme`;
- `filter_evals_to_build` -- drop rows a later build invalidated, per-change rather than by
  version cutoff (`BUILD_CHANGES`, `R/optimizer.R:45`).

A configuration evaluated only in some *other* domain therefore counts as untried here and
will be run again to build in-domain evidence.

### 4.2 The transform and the weights

```
z_ij   = atanh( clamp(score_ij, +/-0.999) )            Fisher z        (.to_z)
1/w_ij = 1 / max(n_test_ij - 3, 1)  +  1 / cap         weight          (.eval_weight)
```

`i` indexes a configuration, `j` a trial; `cap` is `settings$weight_cap`. The clamp keeps a
perfect correlation finite. **Both the estimator below and the surrogate of §8 use these, from
this one helper** — see `LESSONS.md` #35.

The weight is two precisions in series, because two independent errors corrupt a score and they
behave differently as the trial grows:

| | variance | as `n_test` grows |
|---|---|---|
| sampling error of the correlation | `1/(n_test - 3)`, exactly, on the z scale | shrinks |
| config × trial interaction | `1/cap` | **does not shrink** |

The second term is why the textbook `n_test - 3` alone is wrong here. A bigger trial measures
*that trial's* idiosyncrasy more precisely; it does not make the idiosyncrasy smaller. Without
the cap the weight spans 7 to 1434 within a single run (205×) and the largest trials decide
every estimate; at `cap = 50` the span is 8×.

`cap` is the weight an unboundedly large trial would receive, so it reads directly as a trial
size — **50 means "no trial counts for more than one with 53 usable accessions"**. `Inf`
restores the uncapped textbook weight.

> **How `cap` is identified, and what it cost.** `lmer` reports `sd_resid` at unit weight, so if
> `n_test - 3` were the whole variance story the fitted `sd_resid` would be `1`. Across the seven
> runs below it spans **0.096 to 4.72**. Solving `sd_resid = 1` per run
> (`../dev/calibrate_weights.R`) gives `C ≈ 33, 41` for the CV0 runs and `C ≈ 608, 1058` for the
> CV00 runs; three runs cannot inform it at all because their `n_test` is near-constant, which
> makes the fit invariant to the cap. One conservative 50 was chosen over two per-scheme
> numbers. The cost is in §4.6 and `LESSONS.md` #36: at 50, `sd_config` fits to **0** in both
> 24-trial runs. `sd_config = 0` collapses every configuration's BLUP onto the grand mean, so
> that fit supplies nothing to rank on — §4.6. That is a property of the fit on those rows, not
> a finding that the configurations are alike: on `80cfbd7c27e9` the zero is singleton dilution
> and replication recovers a non-zero component (§11), on `e7517df121d5` it does not.

### 4.3 The fit

```r
lme4::lmer(.z ~ 1 + (1 | trial_id) + (1 | config_hash), data = d, weights = .w)
```
(`.blup_scores`, `R/optimizer.R:257`.) The implied model is

```
z_ij = mu + c_i + t_j + e_ij
```

- `mu` -- the grand mean, the typical z score across everything in the slice;
- `c_i` -- **the configuration effect**, the quantity being optimized. `sd_config` is its SD;
- `t_j` -- **the trial effect**, how easy or hard trial `j` is. `sd_trial` is its SD;
- `e_ij` -- everything left over. `sd_resid` is its SD *at unit weight*.

Two things this model does **not** have, both deliberate and both consequential:

- **No `config x trial` interaction term.** It cannot be fitted without evaluating the same
  pair twice, which the design forbids (a duplicate pair is counted as a defect by
  `duplicate_cells`). So "this configuration suits this kind of trial" is confounded into
  `e_ij`, and `sd_resid` is an upper bound on measurement noise, not an estimate of it.
- **No fixed effects for trial descriptors.** Program, year, location and test-set size explain
  essentially none of `t_j` (`LESSONS.md` #19), so the trial is blocked on its own identity.

The fit is allowed to warn. `lmer` warns routinely here because `sd_config` sits near the
boundary against `sd_trial`, and discarding the fit for a warning once silently disabled an
entire scheduling tier (`LESSONS.md` #30). Warnings are muffled and recorded in
`estimator_note`; only an error gives up, and then `mean_score` falls back to the pooled
estimate and `se` is `NA` for every configuration.

### 4.4 What comes out

- `mean_score_i = tanh(mu_hat + c_hat_i)` -- the BLUP, on the correlation scale. **This is the
  number the search ranks on**, everywhere. It is a posterior mean, so one configuration's
  `mean_score` exceeding another's says that the posteriors are centred differently, nothing
  about whether the data can tell them apart; §4.5 says how much separation that ordering is
  worth.
- `se_i` -- the posterior SD of `c_hat_i`, on the z scale.
- `z_mean_i = mu_hat + c_hat_i` -- the same estimate as `mean_score`, before the `tanh`, so that
  it and `se_i` are on one scale. The in-play rule of §7 forms its interval here and transforms
  the endpoints; adding a z-scale `se` to a `tanh`'d mean biased the bound in favour of whatever
  was already winning (`LESSONS.md` #43).
- `pooled_i` -- the weighted mean, kept alongside so the two can be compared.
- `n`, `n_ok` -- rows and usable rows.

Of these, **`se` is read by the in-play rule of §7 and by nothing else directly** -- but that
rule now governs two decisions, not one: which configurations earn further evaluations, and which
are parents for crossover and mutation (§8.5). Until 2026-11-04 parenthood ranked on `mean_score`
alone, so a configuration could be judged good enough to deepen and not good enough to breed
from, with nothing giving a reason for the asymmetry (`LESSONS.md` #42). The EI reference point
(§8.4) and the reported incumbent still rank on `mean_score` alone. §7 is where the framing of §2
does its concrete work.

### 4.5 Shrinkage, and the precision replication can buy

The BLUP is a compromise between a configuration's own data and the prior that it is average.
Approximately (exact only for a balanced design; the code says so at `R/optimizer.R:266`):

```
c_hat_i  =  lambda_i * ( zbar_i - mu_hat )

                    W_i / sd_resid^2
lambda_i = ---------------------------------      W_i = sum of that config's weights
           1/sd_config^2 + W_i / sd_resid^2

se_i     ~  sqrt( 1 / ( 1/sd_config^2 + W_i/sd_resid^2 ) )   <=  sd_config
```

Read off three facts:

- **`se` can never exceed `sd_config`.** With no data at all, a configuration's uncertainty is
  the prior, and the prior is `sd_config`.
- **The data half-outweighs the prior at**

  ```
  n* = ( sd_resid / sqrt(w_bar) )^2 / sd_config^2          w_bar = median weight
  ```

  evaluations -- at `n = n*`, `lambda = 1/2` and `se = sd_config/sqrt(2)`. `n*` is the exchange
  rate between replication and precision on this target, because

  ```
  se_i / sd_config  =  1 / sqrt( 1 + n_i/n* )
  ```

  and a configuration can be evaluated at most once per trial, so `n_i <= N`, the number of
  usable trials in the universe. The **best precision the target allows**, at any budget, is
  therefore `se_min/sd_config = 1/sqrt(1 + N/n*)`.

  **This is a bound on precision, not a verdict on whether the question can be answered.**
  Configuration effects are always estimated with uncertainty here; `n*` says how much of that
  uncertainty a full budget can remove, and `n*/N` is what decides it. At `n* = 3` against
  `N = 24`, `se` can be driven to 0.33 of the prior, so a configuration a full `sd_config` above
  the mean would sit three posterior SDs clear of it. At `n* = 25` against `N = 24`, `se`
  reaches only 0.71 of the prior, so even a saturated budget leaves most configurations'
  posteriors overlapping almost entirely: the search can still name the configuration that currently looks
  best, and that is still the right one to exploit, but it is a guess the data cannot
  distinguish from many alternatives. Both regimes occur in this project (§4.6), and the
  practical difference between them is how much weight to put on the identity of the leader
  versus the behaviour of the whole contender set.
- **A configuration's rank is not its raw mean.** At low replication the BLUP is nearly the
  grand mean, which is the intended protection against a lucky single trial -- but it also
  means that ranking at low replication is ranking small residuals of a shrunken quantity.

### 4.6 Measured values, all seven runs

Final state of each run, re-fitted with the optimizer's own functions at the live
`weight_cap = 50` (`../dev/measure_churn.R`, `../dev/calibrate_weights.R`). The `sd_config` the
*uncapped* weight produced is shown beside it, because that difference is the whole of §4.2.

Every `sd_config` in this table is the **fitted** value -- the REML estimate on everything the
run evaluated -- which is the number the shrinkage used through 2026-11-03 and the number §11
argues is the wrong population. The reference-population estimate beside it, and what it does to
the `in play` column, is in §11's table; the two differ by a factor of 2 to 4 where both exist
and by everything where the fitted one is 0.

| run | scheme | trials usable | `sd_trial` | `sd_config` | (uncapped) | `sd_resid`/sqrt(w) | `n*` | in play |
|---|---|---|---|---|---|---|---|---|
| `80cfbd7c27e9` | CV00 | 24 | 0.101 | **0.0000** | 0.0135 | 0.086 | -- | 2147 / 2147 |
| `46b0d4b39c79` | CV00 | 4 | 0.178 | 0.0200 | 0.0200 | 0.020 | 1.0 | 485 / 660 |
| `f92e28ec6415` | CV00 | 4 | 0.182 | **0.0000** | 0.0000 | 0.009 | 4e10 | 7343 / 7343 |
| `7d75f214584a` | CV00 | 9 | 0.061 | 0.0285 | 0.0475 | 0.079 | 7.7 | 252 / 1147 |
| `34e15bd03838` | CV0 | 9 | 0.137 | 0.1001 | 0.1727 | 0.178 | 3.2 | 105 / 350 |
| `e7517df121d5` | CV0 | 24 | 0.223 | **0.0000** | 0.0305 | 0.167 | -- | 705 / 705 |
| `9a2664962c01` | CV0 | 8 | 0.234 | 0.0078 | 0.0078 | 0.013 | 2.9 | 429 / 1023 |

**`sd_config = 0` is a distinct regime, not a large `n*`.** `n*` is a division by `sd_config²`,
so it is undefined there. What the fit is reporting is not "this much more replication would
help" but "on these rows I can see no configuration differences at all": every configuration's
BLUP shrinks to exactly the grand mean, so all of them tie at the maximum, **all** of them
satisfy the contender test, and the caps then hand replication and parenthood an arbitrary few.
Read it as a statement about this fit on this population rather than about the configurations --
§11 separates the two readings by re-fitting on better-replicated subsets, and the reference
`sd_config` it describes is what rescues exactly this case: on the three runs where the fitted
value is 0, the reference estimate supplies the standard errors the rule needs.

**Why two runs have only 4 usable trials.** `46b0d4b39c79` and `f92e28ec6415` are CV00 runs on a
4-trial domain, and that domain was **forced, not chosen**: under CV00 every other trial in the
target was infeasible, and those four were the only ones that were not. CV00 excludes any
accession appearing in the focal trial's genotype file, which on most of these trials leaves too
few training accessions to build a GRM at all (`scripts/CLAUDE.md`, "Platform mismatch"). So the
small domain is a property of the *data under that masking scheme*, not a configuration mistake.
The two runs are therefore not comparable with the 9- and 24-trial runs on anything that scales
with the number of trials: §4.5's precision bound `1/sqrt(1 + N/n*)` has `N = 4` for them, and
`f92e28ec6415`'s `n* = 4e10` is not a finite replication requirement but the signature of
`sd_config` sitting at the boundary. Read them as the hardest available case, not a degenerate
one.

Three of seven runs are in the `sd_config = 0` regime at `cap = 50`; one (`f92e28ec6415`) was
there uncapped too, and two (`80cfbd7c27e9`, `e7517df121d5`) are put there by the cap. `80cfbd7c27e9` keeps a
non-zero `sd_config = 0.0105` at its own calibrated cap of 1058, so for that run the collapse is
the cap's doing and not the data's — see `LESSONS.md` #36 for why one conservative number was
chosen anyway.

**Under the calibrated weight the field is less resolved, not more.** On exactly the same rows
and the same pools, the *uncapped* in-play counts are 8 (`80cfbd`), 5 (`e7517`), 63 (`7d75`),
45 (`34e15`) and 301 (`9a26`); at `cap = 50` they are 2147, 705, 252, 105 and 429. Capping
shrinks the configuration spread the fit can see faster than it shrinks `se`, so more
configurations' uncertainty still covers the leader. The reading is not that the cap lost
information but that the uncapped weighting had been making the field look better resolved than
these data support, by letting the largest trials set nearly the whole estimate (§4.2). Which of
the two is right is the calibration question of §4.2, and uncapped `sd_resid` exceeded 1 -- rows
carrying variance the weight does not model -- in four of the seven runs (4.72, 2.59, 1.32,
1.11).

And `se` as a function of replication, as a fraction of `sd_config` (median per bin), for the
runs where `sd_config > 0`:

| run | `n*` | n=1 | n=2 | n=3 | n=4-5 | n=6-8 | n=9-12 |
|---|---|---|---|---|---|---|---|
| `7d75f214584a` | 7.7 | 0.94 | 0.90 | 0.86 | 0.82 | 0.75 | 0.69 |
| `34e15bd03838` | 3.2 | 0.87 | 0.80 | 0.73 | 0.67 | 0.58 | 0.52 |
| `9a2664962c01` | 2.9 | 0.86 | 0.76 | 0.69 | 0.63 | 0.56 | -- |

The `n*` column still predicts the row: the larger it is, the slower `se` falls away from the
prior, and `1/sqrt(1 + n/n*)` tracks each line to within a few percent.

`sd_trial` exceeds `sd_config` in **every** run under the calibrated weighting, by factors from
1.4 to infinity — it was six of seven uncapped. That is the central fact about this objective and
the reason for every mechanism in §5 and §6.

---

## 5. Choosing the next configuration

`choose_config` (`R/optimizer.R:511`) returns the first of these that applies.

### Phase 1 -- seeds

The five Predictathon submissions, scheme-specific, each run until it has a row in this run's
slice. Worker `w` takes the `((w-1) mod length(undone)) + 1`-th un-done seed, so six workers
start six different seeds instead of all six starting Prediction1 (selection sees only what has
*finished*, never what is in flight).

### Phase 2 -- replication

`.replication_backlog` (`R/optimizer.R:384`) builds two tiers from the run's slice:

- First, drop any configuration already run on **every** trial in the universe: nothing is left
  to learn from it.
- **`base`** -- every configuration with fewer than `config_replication` (2) **rows**. Failures
  count, so the tier always drains.
- **`extra`** -- configurations at or above the floor **that are currently contenders** (§7).

The caller then rations:

```r
pool <- if (length(bl$base)) bl$base
        else if ((nrow(evals) + w) %% every == 0L) bl$extra
        else character()
```

- **`base` preempts `extra` entirely.** While any configuration is below the floor, no
  contender earns a deepening evaluation. Every newly proposed configuration enters at 1 row,
  so with several workers proposing, `base` is refilled continuously.
- `extra` is rationed to one worker in `replicate_every` (3), keyed on the store's row count so
  a lone worker alternates over time and N workers split without lockstep.
- Within whichever pool is chosen, worker `w` takes entry `((w-1) mod |pool|) + 1`.

### Phase 3 -- random initialization

Until `n_random_init` (25) distinct configurations have a finite `mean_score`, draw a fresh
random configuration whose hash is not already in the slice. Logged as `source = "random_init"`.

### Phase 3b -- the persistent reference rate

One step in `reference_every` (20) draws a fresh random configuration **however mature the run
is**, logged as `source = "random_reference"`. Keyed on the store's row count plus the worker
index, exactly as `replicate_every` is, so one worker alternates over time and N workers split
at any instant without lockstep.

This exists because `n_random_init` is a one-time gate on a cumulative count: once it is passed
the reference population never grows again, and the `sd_config` estimated from it (§11) is frozen
at whatever the first 25 configurations could support -- which on one real run was not enough to
support anything. A rate, not a budget.

It sits after replication so it cannot preempt the `config_replication` floor, and before
acquisition so the 5% it costs is taken out of exploitation. That is the honest accounting: its
purpose is to buy an unbiased view of the space, and exploitation is what pays for it.

### Phase 4 -- surrogate-guided

Fit a surrogate, propose candidates, pick the one with the highest EI. Detailed in §8.

### Phase 5 -- fallback

If the surrogate cannot be fitted at all, a fresh random configuration.

> **Where the budget actually went.** Evaluations per distinct configuration, measured:
> 1.7 (`7d75`), 3.2 (`80cfbd`), 3.4 (`f92e`), 4.2 (`e7517`), 4.3 (`34e15`), 5.2 (`46b0`),
> 6.0 (`9a26`). Configurations that reached **every trial that was usable in their run**:
> 10 of 2147 (`80cfbd`), 1 of 705 (`e7517`), 7 of 1023 (`9a26`), 15 of 350 (`34e15`),
> 10 of 1147 (`7d75`). And the share of usable evaluations spent on configurations that never
> got past 2 of them -- the unrationed `base` floor -- was 60% in `7d75` and 61% in `46b0`,
> against 13-15% in `9a26` and `e7517`.

---

## 6. Choosing the trial

`choose_trial` (`R/optimizer.R:451`), given the configuration just chosen.

1. **Exclude trials this configuration has already been run on**, plus those other workers are
   running it on right now (`active_claims`). The pipeline is deterministic in
   `(config, trial, scheme)`, so repeating a pair recomputes a known score and double-counts it.
2. **Compute the trial replication target**: `.trial_target` = `trial_replication` +
   `floor(sqrt(n_scored)/5)`, where `n_scored` is the number of distinct configurations with a
   finite score in the slice. The square root is deliberate: under a linear schedule the number
   of *new* trials per unit budget flattens to a constant and the search stops seeing new
   environments.
3. **Serve the trial backlog first**: trials that have been started but carry fewer than that
   many distinct configurations. Worker `w` takes entry `((w-1) mod |backlog|) + 1`. This is
   the precondition for `trial_id` being usable as a surrogate feature at all -- at one
   observation per trial an `rpart` split on `trial_id` memorizes noise to an in-sample R^2 of
   0.98 (`trial_feature_usable`, `R/surrogate.R:111`).
4. **Otherwise draw a fresh unseen trial** from the eligible catalogue, up to 5 tries.
5. If every draw was already seen, raise `trials_exhausted`. `optimizer_step` then re-runs
   `choose_config` with `replicate = FALSE` -- spending the step exploring rather than idling,
   and deliberately not counting it toward `max_sample_fail`, so a saturated contender set
   cannot halt the run.

A trial in the universe that `sample_real_trial` always rejects -- too few genotyped
accessions, or no focal-trait observations under the trait-*name* filter -- is skipped silently
on every draw. The consequence is not only wasted draws: because the universe size is used as
the "domain-covered" test in §5, **an unreachable trial means no configuration is ever dropped
from the backlog for having covered the domain.**

> **Measured.** Run `9a2664962c01` pinned a 15-trial universe; only 14 trials were ever
> evaluated. The missing one is trial 10154 -- the same trial named in the code comment at
> `R/data_access.R:535` for having been rejected on every draw for 2,471 evaluations.

---

## 7. What counts as still in play: one rule, two uses

```r
.in_play <- function(agg, z = 0.25, delta = 0, rule = "racing") {
  a <- agg[is.finite(agg$mean_score) & is.finite(agg$se %||% NA_real_), , drop = FALSE]
  if (!nrow(a)) return(character())
  zm <- if ("z_mean" %in% names(a)) as.numeric(a$z_mean) else rep(NA_real_, nrow(a))
  zm[!is.finite(zm)] <- .to_z(a$mean_score[!is.finite(zm)])   # exact fallback
  ucb <- tanh(zm + z * a$se)
  bar <- if (identical(rule, "point")) max(tanh(zm)) else max(tanh(zm - z * a$se))
  keep <- ucb + delta >= bar
  a$config_hash[keep][order(ucb[keep], decreasing = TRUE)]
}
```
(`R/optimizer.R`, with `.contenders()` taking the first `k` of it.)

In words, in order:

1. **Drop every configuration without an `se`.** `se` comes only from the random-effects fit, so
   if that fit does not run there are no contenders, `extra` is permanently empty, and
   `replicate_every` has nothing to ration (`LESSONS.md` #30). `get_parents()` carries a floor
   for exactly this case, so the parent pool does not silently empty with it (#42).
2. **Form the claim**, `ucb_i = mean_i + contender_z * se_i`: an optimistic view of
   configuration *i*.
3. **Form the bar.** Under `contender_rule = "racing"` it is `max_j(mean_j - z*se_j)`, the best
   value that can be *defended* about any configuration. Under `"point"` it is `max_j(mean_j)`.
4. **In play** if `ucb_i + resolution_delta >= bar`.
5. **Order by UCB, descending.** The first 8 earn contender replication (§5 Phase 2); the first
   `n_parents` = 32 are the parents of crossover and mutation (§8.5).

### Why the bar is a lower bound and not a point estimate

`max_j(mean_j)` is the maximum of a set of noisy estimates — the most upward-biased order
statistic in the set. One thinly replicated configuration with a lucky draw therefore sets the bar
for everyone, **including for a well-replicated leader, which it can evict.** It is also
unprincipled on its own terms: the candidate is charged `z` for its uncertainty and the leader is
given a free pass on its.

The racing bar charges both sides the same `z`, which is the comparison a racing algorithm
eliminates on: a candidate goes out only once its interval stops overlapping the leader's. With
A at 0.50 ± 0.02, B at 0.45 ± 0.08, C at 0.48 ± 0.01 and a one-trial D at 0.54 ± 0.15, at
`z = 0.5`:

| | bar | in play |
|---|---|---|
| point | 0.54 (D's mean) | D alone — **A, the leader, is evicted** |
| racing | 0.49 (A's own lower bound) | D, A, B |

Note which one racing admits: **B, the uncertain one**, not C, which is precisely estimated and
precisely second. It keeps what is unresolved, not what is merely close. And D joins the set
without moving the bar at all, because setting the racing bar takes a high mean *and* a small `se`.

Two consequences of the choice:

- **`z` does double duty under racing** — it inflates the claim and lowers the bar — so the set
  grows about twice as fast in `z`, and `z_racing ~ z_point / 2` gives comparable set sizes. At
  `z = 0` the two rules are identical. This is why the default `contender_z` is 0.25 rather than
  the 1 that the point rule ran with.
- **The comparison happens on one scale.** The estimate and its `se` are both on z, where the
  fit produced them, and only the interval *endpoints* are transformed:
  `ucb_i = tanh(z_mean_i + z*se_i)`, `bar = max_j tanh(z_mean_j - z*se_j)`. Until 2026-11-04 the
  code added a z-scale `se` to a `tanh`'d mean, which is on neither scale and overstates the
  bound by up to a quarter at the top of the range (`LESSONS.md` #43). Because `tanh` is
  monotone, the fix changes no in/out decision for `z` alone -- but `delta` is a difference in
  predictive *accuracy*, which only means anything on the correlation scale, so the comparison
  has to be there for `delta` to be added to it at all. `z_mean` comes from the random-effects
  fit; `atanh(mean_score)` is the exact fallback for a caller that supplies only `mean_score`.

### `resolution_delta`: the one term the search cannot shrink

Everything else in the test is denominated in `se`, and `se <= sd_config`, which the search itself
drives down (§11). `resolution_delta` is the difference in predictive accuracy actually worth
resolving, chosen from *outside* the run — so a configuration keeps its place when it misses the
bar by less than that, and membership flips only on evidence that moves further than the
resolution anyone cares about. A deadband, in other words, and the only term here that does not
move when the population homogenizes.

> **Measured** (`dev/sweep_contender_rule.R`, racing at `z = 0.25` with the reference
> `sd_config`, pooled over the runs whose fitted `sd_config` is non-zero). `jump` is the median
> |log₂| change in set size between consecutive checkpoints: 0 means a set that never moves, 1
> that it typically doubles or halves. `usable` is the share of checkpoints where the set is at
> once big enough to fill the contender cap of 8 and no more than a tenth of what has been scored.
>
> | `resolution_delta` | in play (median) | % of scored | `jump` | `usable` |
> |---|---|---|---|---|
> | 0.000 | 6.0 | 1.9 | 0.70 | 0.23 |
> | **0.002** | **11.5** | **5.1** | **0.63** | **0.23** |
> | 0.005 | 32.0 | 11.5 | 0.33 | 0.19 |
> | 0.010 | 79.5 | 23.3 | 0.28 | 0.07 |
> | 0.020 | 116.5 | 59.5 | 0.24 | 0.04 |
>
> **The default is chosen for set size, not for stability.** 0.002 roughly doubles the median
> set, from just under the cap of 8 to just over it, and buys only 9% off `jump` while doing it.
> The real stability is at 0.005 and above — but a set of 80 or 117 can hardly move *because* it
> holds a quarter to three fifths of the field, which is the sense in which `jump` rewards
> vacuity and must always be read beside the set size. `usable` says the same thing in one
> number, and it is flat from 0 to 0.002 and falls thereafter.
>
> Per-run, two of the seven end up with a median set **below** the cap of 8 (`f92e28ec6415` at
> 1.5, `e7517df121d5` at 2.5): there the field is resolved at this resolution — honestly on the
> first, on a prior the next section says to distrust on the second — the `extra` replication
> tier is small by design, and the `n_elites` floor of §8.5 is what keeps a parent pool in
> existence.

Read it against `sd_config`, which the reference population puts at 0.013 to 0.056 across these
runs. **As `delta` approaches `sd_config` everything is in play, and that is the honest reading
rather than a failure:** no configuration then differs from another by an amount worth having, and
the right response is to stop searching, not to search harder.

### What happens when a configuration falls out of the top 8

This is worth stating plainly because the consequence is near-total and silent.

- It is not in `extra` (that tier is `n_eval >= floor` **and** in the first 8 of the in-play set),
  and not in `base` (it is above the floor). **Its priority for further evaluation becomes zero.**
- It can still be a **parent**, if it is within the first `n_parents` = 32. That is the one thing
  the §42 unification softened: under the old design, parents were the mean-ranked top 8, so
  leaving *that* set cost parenthood as well, and the two sets were not even the same 8.
- Acquisition will not pick it up: `propose_candidates` prefers candidates whose hash is not
  already in the slice, and `fresh_random` avoids them too. The one exception is when fewer than 5
  novel candidates survive de-duplication, in which case already-evaluated configurations re-enter
  the pool.
- **Nothing is remembered.** There is no record that it was ever a contender, no partial credit,
  no accumulated priority. The set is recomputed from scratch from the current fit at every step,
  so a configuration re-enters only by having a top-8 UCB again under a later fit. This is where
  the analogy with racing algorithms stops: there, elimination is final and that finality is what
  buys the sample-complexity guarantee. Here nothing is eliminated and nothing is settled, so only
  the comparison transfers, not the guarantee.
- It is frozen at whatever replication it had. No new evidence about it is ever gathered; its BLUP
  still moves, but only because other configurations' data moved the grand mean and the variance
  components. It gets a second chance only on someone else's data, never on its own.

The cap interacts badly with a wide in-play set. Being in play means sitting within one
`contender_z * se` of the bar — so when 301 of 1023 configurations are in play (run
`9a2664962c01` under the old settings), the eight that get deepening evaluations are the top 8 of
301 that the fit itself says it cannot tell apart. That ordering is re-derived after every refit.
Each of the eight gets one extra evaluation and is replaced.

> **Measured (the treadmill).** In `9a2664962c01`, replayed at 14 checkpoints about 464
> evaluations apart, the leader at one checkpoint sits at median rank **19** at the next
> (quartiles 9-79; observed ranks 1, 1, 3, 9, 10, 18, 19, 34, 36, 79, 88, 150, 496, out of
> 32-1023 scored configurations). It does not merely lose the lead -- it leaves the top 8, and
> with it all priority for further evaluation.

The mechanism behind that is the shrinkage of §4.5, not bad luck. The only configurations
receiving extra evaluations are the ones currently estimated highest, which is partly why they are
there; an additional evaluation drawn from a population whose `sd_config` is near zero moves the
BLUP back toward the grand mean. **Replication, as scheduled, systematically removed
configurations from the contender set rather than accumulating evidence on any one of them** — and
the reference `sd_config` of §11 is the change that addresses the cause rather than the symptom.

---

## 8. The surrogate and the acquisition function

### 8.1 The rows the surrogate sees

`choose_config` runs the evaluation slice through **`.fisher_weights()`, the same helper the
estimator uses**, keeps the rows it marks `.usable`, and then **drops every row whose trial
carries fewer than 2 distinct configurations** -- the same safety condition as
`trial_feature_usable`.

**The surrogate uses the statistics of §4.** One usability floor (`min_test_acc`), one weight
(`.eval_weight`, capped at `weight_cap`), one scale (Fisher z), applied by one function, to both
the aggregate and the forest. Concretely:

- the response is `.z`, not `score` -- so the additive trial effect that MERF and BLOCKED
  subtract is additive on a scale where that is coherent (a `+0.3` trial shift applied to a
  configuration at `0.85` leaves the correlation range; on z nothing does);
- `.w` is passed as `rpart` **case weights**, so an anova split minimises *weighted* SSE and a
  node's value is the *weighted* mean, and to the `lmer` that estimates the trial effect;
- `best` in the EI comparison is `atanh(max(mean_score))`, recovering the BLUP the fit produced,
  so incumbent and candidates are finally on one scale.

Note that EI is **not** invariant under `atanh`: it integrates the *size* of the improvement, so
on z a gain at the top of the range counts for more than the same gain in the middle. That is
the intended reading, and `ei_xi` is now a z-scale margin.

Until 2026-11-03 none of this was true: the filter was `is.finite(score)`, the response was the
raw `score` column, and there were no weights anywhere in the surrogate -- so a correlation over
12 accessions and one over 900 were the same observation to the forest while `aggregate_scores`
weighted them 9 to 897. See `LESSONS.md` #35 and #36.

### 8.2 Three modes, most structured first

- **MERF** (default, `surrogate_method = "merf"`): one row per `(config, trial)`. The trial
  effect is a shrunken random effect fitted *outside* the forest: fit `lmer(y ~ 1 + (1|trial))`,
  fit the forest to `y - b`, re-estimate `b` from `y - f(X)`, iterate 3 times
  (`fit_surrogate_merf`, `R/surrogate.R:66`). A many-level grouping factor then costs one
  variance parameter instead of competing for tree structure, and the forest's prediction is
  already trial-marginal.
- **BLOCKED**: the same rows with `trial_id` as an ordinary forest feature, averaged over all
  observed trials at prediction time (`predict_surrogate_marginal`).
- **POOLED**: one row per configuration mean, when the design cannot support a trial effect.

Each falls through to the next if it cannot be fitted, and if none can, the step spends itself
on a fresh random draw. The forest itself is 300 `rpart` trees, each on a row bootstrap and a
random `sqrt(p)` subset of features.

A mode needs `settings$surrogate_min_obs` (25) rows to be attempted at all. **That counts the
rows of whichever mode is being tried, and the modes count differently**: MERF and BLOCKED get
one row per `(config, trial)`, POOLED gets one row per configuration. So a threshold above
`n_random_init` leaves POOLED permanently unfittable, since Phase 3 guarantees only
`n_random_init` scored configurations before acquisition starts. Until 2026-11-03 this threshold
*was* `n_random_init`, which is why the two defaults are equal and why the coupling went
unnoticed: it silently guaranteed POOLED could always fit (`LESSONS.md` #39).

### 8.3 Prediction, and the uncertainty floor

`predict_surrogate` returns the mean across trees and the SD across trees, then:

```r
sd_floor <- 0.25 * stats::sd(model$y_obs)
sd_hat   <- pmax(sd_hat, sd_floor)
```

The stated intent is that "a candidate the forest happens to agree on is not treated as known
with certainty". The effect depends entirely on how many candidates the floor binds for:

- If it binds for **some**, it raises those candidates' EI -- exploration, as intended.
- If it binds for **all**, every candidate has the *same* `sd_hat`, and then EI is a strictly
  increasing function of `mean_hat` (`dEI/dmean = Phi(z) > 0`). **The acquisition function
  becomes greedy exploitation of the surrogate mean, with no exploration term whatsoever.**

This is one of two routes to the same degeneracy; §8.4 is the other, and §8.4's repair does not
touch this one.

### 8.4 The reference point

```r
.marginal <- function(f) if (mode == "blocked") predict_surrogate_marginal(surro, f, trials)
                         else                   predict_surrogate(surro, f)
pr   <- .marginal(cfeats)                       # the candidates
best <- max(.marginal(ifeats)$mean)             # the configurations already evaluated
ei   <- expected_improvement(pr$mean, pr$sd, best, settings$ei_xi)
```

**Both sides of the comparison come from one estimator.** `best` is the forest's own
trial-marginal prediction at the features of the configurations already evaluated -- not
`max(mean_score)`.

It was `max(mean_score)` until 2026-11-03, and that is a scale error with a specific
consequence. A BLUP is shrunk toward the grand mean in proportion to `sd_config` (§4.5); a
forest prediction is not shrunk at all. So `mean_hat - best` compared two different estimation
scales, and the gap grew as `sd_config` fell. Writing EI out,

```
EI = (m - best - xi)·Phi(z) + sd·phi(z),      z = (m - best - xi)/sd
```

driving `best` down sends `z -> +inf`, so `Phi(z) -> 1`, `phi(z) -> 0`, and `EI -> m - best - xi`.
The ranking does not flatten -- it collapses onto `m`. **`sd·phi(z)` is the only term that ever
lets EI prefer a candidate the forest is unsure about, so losing it makes acquisition plain
greedy exploitation.** In the limiting case it is stark: when `sd_config` fits to 0, every BLUP
*is* the grand mean, so `best` was literally the grand mean of the whole run.

The repair cannot run the other way. Shrinking `m` to the posterior scale needs the candidate's
own shrinkage factor `lambda_i = W_i/sd_resid² / (1/sd_config² + W_i/sd_resid²)`, and an
unevaluated candidate has `W_i = 0`, hence `lambda = 0` and every candidate collapses to
`mu_hat`. The posterior scale has no room for a configuration with no data, so harmonisation has
to go toward the unshrunken scale.

> **Measured, before and after**, re-fitting the real surrogate on each run's final rows at
> `weight_cap = 50` and proposing a real candidate pool:
>
> | run | candidates | at the sd floor | `best` old | `best` new | Spearman(EI, mean) old | new |
> |---|---|---|---|---|---|---|
> | `9a2664962c01` | 129 | **129 (100%)** | 0.3525 | 0.3472 | **1.000** | **1.000** |
> | `80cfbd7c27e9` | 172 | 102 (59%) | 0.1211 | 0.2239 | 0.682 | 0.260 |
> | `e7517df121d5` | 173 | 89 (51%) | 0.3516 | 0.5228 | 0.766 | -0.017 |
> | `34e15bd03838` | 170 | 12 (7%) | 0.3671 | 0.3479 | -0.027 | 0.075 |
>
> The size of the correction tracks the shrinkage exactly as the algebra says it should. On the
> two runs where `sd_config` fits to 0, `best` was the grand mean and rises by 0.10 and 0.17,
> and EI stops tracking the predicted mean (Spearman 0.68 -> 0.26 and 0.77 -> -0.02). On
> `34e15bd03838`, where `sd_config = 0.1001` and shrinkage is mild, `best` moves by -0.02 and
> little changes. **On `9a2664962c01` nothing changes at all**, because the sd floor binds for
> 129 of 129 candidates and under a uniform `sd` EI is monotone in `m` whatever `best` is --
> §8.3's route to the same place, which this repair does not address.
>
> `best` is a prediction at a configuration's own features, so it is partly in-sample. Measured
> against the unshrunk pooled means, the forest smooths far harder than it memorises:
> `max(prediction)` is 0.22-0.52 where `max(pooled z)` is 0.70-0.93, and restricting to
> configurations with `n_ok >= 3` barely moves it (0.4948 against 0.5228 on `e7517df121d5`). An
> out-of-bag reference would remove the residual bias; whether it is worth the bookkeeping is
> open -- `LESSONS.md` #37.
>
> In every run the forest predicts systematically lower for fresh random draws than for elite
> recombinations (medians 0.429 vs 0.459/0.474 in `9a2664962c01`; 0.118 vs 0.129/0.131 in
> `80cfbd7c27e9`), because the elite neighbourhood is where the high observed scores are. Under
> greedy selection a random draw therefore **cannot** win, and the only remaining source of
> novelty is mutation. That is what the restored `sd·phi(z)` term exists to counteract.

### 8.5 Candidate generation

`propose_candidates` (`R/optimizer.R`) builds 60 crossovers of parent pairs, 60 mutations of
single parents (resampling 1 or 2 of the 6 blocks), and 60 fresh random configurations, then
de-duplicates by hash and prefers the novel ones.

**The parents are the in-play set of §7**, via `get_parents()`: ordered by UCB, capped at
`n_parents` = 32 and floored at `n_elites` = 8. Both bounds are load-bearing. The cap, because
when the estimator cannot separate anything the in-play set is *every* scored configuration, and
unbounded parenthood is then recombination among all of them. The floor, because the in-play set
is *empty* whenever the random-effects fit did not run -- no fit means no `se`, and §7 drops every
row without one -- and an empty parent pool would turn this function into pure random draws with
nothing saying so.

**Crossover takes whole subtask blocks from one parent or the other.** It cannot create a value
neither parent has. So once the parents agree on a block, crossover is a no-op for that block
forever, and the only ways out are mutation -- which touches a given block with probability
about 1/4 per mutant -- and the random draws, which §8.4 shows can be systematically rejected.
Widening the pool from 8 to as many as 32 weakens this but does not remove it: the measurement
below was taken with 8, and the parents are still drawn from one neighbourhood.

> **Measured (lineage collapse).** In `9a2664962c01`, **97% of the 1043 distinct configurations
> evaluated share identical values of five continuous parameters** (`kernel.maf` =
> 0.0760086870309897, `kernel.ridge`, `kernel.max_missing`, `pheno_prep.z_thr`,
> `pheno_prep.ge_bandwidth`); there are 26 distinct `kernel.maf` values among 1043
> configurations. In the candidate pool rebuilt from that run's elites, the 60 crossover draws
> de-duplicated to **2** distinct novel candidates. Eight successive leaders differed only in
> `train_select.secondary_min` (7/9/10/16/18), `primary_min`, `geno_select.method`,
> `lambda_fixed` and `predict_post`.
>
> This is not inevitable: `80cfbd7c27e9` ended with 2095 distinct `kernel.maf` values among
> 2349 configurations and `e7517df121d5` with 595 among 712.

---

## 9. Churn: a definition

**Churn is turnover among candidate configurations without directional improvement and without
the set of plausible bests ever narrowing.** Three series, measured over a sliding window of
evaluations; all three conditions must hold, because each alone is normal.

| component | metric | normal early | churn |
|---|---|---|---|
| **Turnover** | `1 - Jaccard(top-k now, top-k Delta evaluations ago)`, `k = n_elites` | high | high |
| **Direction** | slope of the leader's `mean_score` against evaluation count, in units of its own `se` | positive | ~0 |
| **Resolution** | the **uncapped** in-play count, `.in_play()` at the live rule and `z` (§7); secondarily the leader's trial coverage and `max(n_ok)` among contenders | falling | flat or rising |

Read *resolution* as narrowing, never as certainty: the series asks whether the set of
configurations whose uncertainty still covers the leader is shrinking, not whether any
configuration has been established as best. Turnover with direction is productive search. Low
turnover with a narrow in-play set is convergence. Churn is high turnover, no direction, and a
field that never narrows.

A fourth, cheaper diagnostic says how much narrowing the target even permits:
**`n*` (§4.5) against `N`, the number of usable trials in the universe.** If `n*` is comparable
to or larger than `N`, the bound of §4.5 holds `se` near `sd_config` however the budget is
spent, and a wide in-play set is then the honest report of what the data support rather than a
scheduling failure.

> **Measured.** Replaying each run's ranking over prefixes of its own evaluation stream, at 25
> evenly spaced checkpoints from evaluation 200 to the end (spacing in evaluations: 305
> `80cfbd`, 115 `e7517`, 83 `7d75`, 54 `34e15`, 251 `9a26`), under the live `weight_cap = 50`.
> Tighter spacing raises Jaccard, so the comparison understates the gap: `9a26` churns most
> despite being sampled no more finely than the rest.
>
> | run | leader changes | distinct leaders | median top-8 Jaccard | in-play at end |
> |---|---|---|---|---|
> | `80cfbd7c27e9` | 7 / 25 | 4 | 0.78 | **2147 / 2147** |
> | `e7517df121d5` | 12 / 25 | 9 | 0.60 | **705 / 705** |
> | `7d75f214584a` | 12 / 25 | 9 | 0.60 | 252 / 1147 |
> | `34e15bd03838` | 14 / 25 | 12 | 0.78 | 105 / 350 |
> | `9a2664962c01` | 21 / 25 | 20 | 0.39 | 429 / 1023 |
>
> The two bolded in-play figures are the `sd_config = 0` regime of §4.6, not a resolution
> failure of the kind this section is about: there the BLUPs are all equal, so the whole
> population ties and `in_play` is trivially everything. Read those two rows' turnover columns
> as meaningless rather than good.
>
> In `9a2664962c01` the leader's score rose from 0.3236 to 0.3419 by evaluation 703 and then
> stayed between 0.336 and 0.339 for the remaining 5533 evaluations, while the leader's identity
> changed at 19 of 25 checkpoints and the leader's usable replication was frequently 1-3.
> Turnover high, direction zero, resolution never achieved: churn, by the definition above.
> `80cfbd7c27e9`, on the same code and the same build, did not churn.

---

## 10. Why it churns: the causal chain

Five mechanisms, each with its code anchor and its measurement. They compound; the first two
set the stage and the last three close a loop.

**1. The objective's noise is structural, not reducible.** `sd_trial` exceeds `sd_config` in
**every** run under the calibrated weighting. The trial effect is removed by blocking, but what
remains (`sd_resid`, which also absorbs the unfittable `config x trial` interaction) is 1.0 to
2.8 times `sd_config` per evaluation on the four runs where `sd_config > 0`, and unbounded on
the three where it fitted to zero. That ratio, squared, *is* `n*`. §4.6.

**2. On some targets the attainable precision is low at any budget.** `se <= sd_config` always,
and replication reduces it only in proportion to `n/n*`, so with `n` capped at the universe size
`N` the best attainable `se/sd_config` is `1/sqrt(1 + N/n*)`. Under the uncapped weight that
bound sat near three-quarters of the prior on the 24-trial targets (`n*` of 24.9 and 26.6
against 24 trials): a saturated budget would still have left the configurations' posteriors
heavily overlapping. Under the calibrated weight those two runs instead fit `sd_config = 0`,
which is the same statement with the prior at zero -- the fit sees no configuration differences
to narrow. §4.5, §4.6. Neither version says the configurations are equally good; both say the
evidence obtainable on that target cannot show which is better. For `80cfbd7c27e9` the zero is
**singleton dilution and therefore recoverable by replication** (§11); for `e7517df121d5` it is
not.

**3. The cap turns deepening into a treadmill.** Only the top 8 by UCB earn extra evaluations;
dropping out costs all replication priority instantly (§7, softened only for parenthood since
2026-11-04); and the extra evaluation itself
pulls a top-ranked configuration back toward the grand mean. Measured: the previous
checkpoint's leader sits at median rank 19 at the next.

**4. Acquisition degenerates to greedy exploitation, and the elites become one lineage.**
Two independent routes removed the `sd·phi(z)` exploration term: strong shrinkage understated
`best` (§8.4), and the uniform uncertainty floor equalised `sd` (§8.3). Either one alone makes
EI a monotone function of the predicted mean. The forest then systematically predicts lower for
random draws, so the pick is always an elite recombination. Crossover among agreeing elites is a
no-op, leaving mutation of 1-2 of 6 blocks as the only novelty. Measured: 100% of candidates at
the floor, Spearman(EI, mean) = 1.000, 97% of configurations sharing five identical continuous
parameters, 60 crossover draws collapsing to 2 distinct candidates. §8.3-8.5.

The `best` route was repaired on 2026-11-03 and the floor route was not, so this mechanism is
now half-open: it still fires in full wherever the floor binds for every candidate, which is
exactly the run it was measured on. The proposal step is not the bottleneck either -- ~15 of the
60 mutants touch any given block, and with 2-4 methods per block roughly 10 candidates per
iteration already carry an off-consensus method. **One is evaluated, chosen greedily**, so
fixation happens at selection, not at generation.

**5. The loop closes: a homogeneous population drives `sd_config` down, which strengthens
shrinkage, which strengthens 3 and 4.** See §11.

**And a data-level amplifier.** A focal trial whose test set cannot reach `min_test_acc` costs
a full evaluation and contributes nothing: the row is stored `too_few_overlap` with `score = NA`,
so it counts toward the replication floor and enters no estimate. Where a whole trial sits below
the gate, the effective `N` is smaller than the pinned universe, so the bound in mechanism 2 is
governed by the smaller number -- and nothing in the report states that `N`. Until 2026-11-02
this was worse: a separate, lower feasibility floor of 5 let `5 <= n_test < 10` rows be stored
`ok` with a real score that every estimate then silently discarded -- 36% of the successful
evaluations in one run (§3). One setting now governs both ends. §3, §6, `LESSONS.md` #34.

---

## 11. Why `sd_config` falls, and whether that is a defect

`sd_config` is the SD of the configuration effects **in the population of configurations
currently in the fit**. It is not a property of the search space; it is a property of what the
search has chosen to evaluate.

> **Measured, run `9a2664962c01`** (checkpoints over its 6236 evaluations):
>
> | evaluations | distinct configs | `sd_config` | SD of BLUPs | median `se` | in-play |
> |---|---|---|---|---|---|
> | 200 | 32 | 0.0263 | 0.0182 | 0.0171 | 7 |
> | 1129 | 197 | 0.0175 | 0.0125 | 0.0102 | 142 |
> | 2057 | 362 | 0.0137 | 0.0099 | 0.0080 | 314 |
> | 3450 | 587 | 0.0107 | 0.0077 | 0.0064 | 532 |
> | 4843 | 810 | 0.0091 | 0.0064 | 0.0056 | 543 |
> | 6236 | 1023 | **0.0078** | 0.0051 | 0.0053 | 301 |

Three candidate causes, separated by measurement:

- **Singleton dilution** -- most configurations have 1-2 rows and carry little information about
  `c_i`. *Small on `9a26`, and the whole story on `80cfbd`.* Re-fitting on better-replicated
  subsets only: `sd_config` rises from 0.0078 to 0.0105 on `9a26` (`n_ok >= 6`) -- real, but
  nothing like the factor of 3.4 this section is explaining. On `80cfbd7c27e9` the same
  restriction takes it from **0.0000 to 0.0068** (`n_ok >= 4`) and **0.0134** (`n_ok >= 6`): its
  whole-population zero is dilution by singletons, not an absence of configuration differences,
  and `config_replication` plus contender replication (§7) are the mechanisms that undo it. The
  contrast matters -- on `e7517df121d5` and `f92e28ec6415` `sd_config` is 0 at **every**
  replication threshold, so there the zero is the answer rather than an artefact.
- **Estimation noise early on** -- at 200 evaluations and 32 configurations the estimate is
  itself unstable (`sd_resid` was 0.48 there against 0.16 at the end). *Contributes to the first
  checkpoint or two, not to the trend.*
- **The population genuinely became homogeneous.** *This is the cause.* Splitting `9a26` in
  half and fitting each half separately on configurations with at least 2 rows *within that
  half*: first half `sd_config` = 0.0109, second half `sd_config` = **0.0000** -- the component
  collapses to the boundary. The configurations evaluated late in that run do not differ.
  In `80cfbd7c27e9`, which did not collapse, the same split gives 0.0152 and 0.0138: flat.

**So: is it inevitable?** A falling `sd_config` is the expected signature of a search that is
concentrating, and on its own it is not a fault. What makes it a defect is three things the code
does with that number:

1. **It is the shrinkage prior.** As `sd_config` falls, `lambda` falls, every BLUP moves toward
   the grand mean, and the ranking turns on ever-smaller differences.
2. **It sets the scale of the contender test.** "Still in play" is measured in units of `se`,
   and `se <= sd_config`. The decision rule is denominated in a unit the search itself shrinks.
3. **It feeds back into the collapse that caused it.** More shrinkage understates `best`, which
   (§8.4) makes EI less selective and more greedy, which proposes more elite recombinations,
   which makes the population more homogeneous, which lowers `sd_config` further.

The distinguishing test between convergence and collapse is whether the leader's score is
rising. In `9a2664962c01` `sd_config` fell by a factor of 3.4 while the leader's score did not
move after evaluation 703. That is collapse, not convergence.

### The repair: estimate it on a population the search did not choose

The quantity the code *wants* is the spread of configuration effects in the configuration
**space**. What it had was the spread in the set the search chose to evaluate. Since 2026-11-03
the store records which phase of `choose_config()` proposed each row (`evals.source`,
`LESSONS.md` #38), so the two can be separated: the Phase 3 random draws and the Phase 5 random
fallback are samples of the space, and everything else is a sample of the search's trajectory.
Seeds are excluded -- five hand-picked submissions are not a sample of anything. `sd_config` is
then estimated on the reference rows alone and the configuration effects are re-shrunk on it
(`settings$sd_config_source = "reference"`, `R/optimizer.R::.reference_sd_config`).

> **Measured** (`dev/sweep_contender_rule.R`, 12 checkpoints per run). `sd_fit` is the REML
> estimate on everything evaluated; `sd_ref` is the same estimator on the reference draws alone.
>
> `in play` is the median over checkpoints, as a share of the scored population, under the
> settings through 2026-11-03 (point bar, `z` = 1, Δ = 0, fitted `sd_config`) and under the
> current ones (racing bar, `z` = 0.25, Δ = 0.002, reference `sd_config`). `jump` is defined in §7.
>
> | run | `sd_fit` first → last | `sd_ref` (median) | in play, old → new | `jump`, old → new |
> |---|---|---|---|---|
> | `34e15bd03838` | 0.0856 → 0.1001 | 0.0560 | 26.8% → 9.5% | 0.44 → 0.16 |
> | `46b0d4b39c79` | 0.0230 → 0.0200 | 0.0273 | 66.4% → 29.5% | 0.45 → 1.00 |
> | `7d75f214584a` | 0.0162 → 0.0285 | *none, at every checkpoint* | 42.4% → 7.5% | 0.47 → 0.63 |
> | `80cfbd7c27e9` | 0.0000 → 0.0000 | 0.0130 | **100% → 0.74%** | 0.25 → 0.29 |
> | `9a2664962c01` | 0.0251 → **0.0078** | **0.0246** | 69.5% → 33.2% | 0.40 → 0.67 |
> | `e7517df121d5` | 0.0000 → 0.0000 | 0.0517 | 78.4% → 0.71% | 1.04 → 1.00 |
> | `f92e28ec6415` | 0.0251 → 0.0000 | 0.0250 | **100% → 0.04%** | 0.24 → 0.00 |

Three things that table says:

- **On `9a2664962c01` the reference estimate is flat where the fitted one collapses** — 0.0240 to
  0.0246 across the run against 0.0251 down to 0.0078. That is the hypothesis of this section
  confirmed directly: the fall was the search narrowing, not the space.
- **It exists where the fitted one is zero.** On the three runs whose fitted `sd_config` reaches
  the boundary, every BLUP ties at the grand mean and every `se` is 0, so §7 is inoperative at
  every `z` — there is nothing for the contender rule to read. The reference estimate supplies
  standard errors there, which is what makes replication and parenthood work at all on those
  runs. `80cfbd7c27e9`'s 0.0130 agrees with three independent routes to the same number (the
  uncapped fit 0.0135, the `n_ok >= 6` restriction 0.0134, a first-30-configurations
  reconstruction 0.0130).
- **`jump` gets worse on four of the seven runs, and that is not a regression.** It is the
  artefact §7 warns about, arriving at full strength: a larger set is mechanically steadier, and
  these sets shrank by up to three orders of magnitude. The old settings were not stable, they
  were **vacuous** — 100% of the population in play on two runs and 66-78% on three more, a set
  that cannot move because it already contains everything. At `e7517df121d5`'s new median of 2.5
  configurations, a `jump` of 1.00 is one member changing. What the new numbers buy is a set
  small enough for the cap of 8 to actually serve; what they expose is that step-to-step
  turnover in a *small* set is now the live problem, which is the "nothing is remembered"
  paragraph of §7 and is not addressed by any of this.
- **It is not uniformly a loosening.** On `34e15bd03838` the fitted value (0.100) is *larger* than
  the reference one (0.056), so that run gets a *tighter* prior. That is the honest counter-case,
  and it is the right behaviour: a run whose evaluated set really is more spread out than the
  space it was drawn from should not have its differences taken at face value.

### Three things that had to be added before the repair was safe

**1. The reference population has to keep growing.** `n_random_init` is a one-time gate on a
cumulative count, so when Phase 3 ends the reference population is frozen at 25 configurations
and its estimate never improves. At `config_replication` = 2 that is roughly two rows each --
thin enough on `7d75f214584a` that the reference fit returns a boundary 0 and declines
altogether. Two settings address it: `reference_every` (§5 Phase 3b) draws a fresh random
configuration one step in 20 forever, and `reference_replication` (4, against 2 for everything
else) gives reference configurations extra evaluations through the *rationed* tier of
`.replication_backlog` -- never the unrationed floor, which would let a reference target preempt
the whole search. Neither helps retroactively: on the runs measured here the population is still
the frozen 25.

**2. The pin must change the prior and nothing else.** The first implementation recomputed the
configuration effects from trial-adjusted deviations with the one-way analytic form, which
ignores the uncertainty in the trial effects and runs a few percent low against lme4's
`postVar`. So switching `sd_config_source` switched the `se` *estimator* as well as the prior,
and the two settings were not comparable. Now, writing
`prec_i(s) = 1/s² + W_i/sd_resid²` for the posterior precision at prior SD `s`, one ratio
`R_i = prec_i(fit)/prec_i(pin)` governs both quantities -- the effect scales by `R_i` and the
`se` by `sqrt(R_i)` -- so lme4's own output is *rescaled* and the result is exactly lme4's at
`pin == fit`.

That route needs a conditioning test, not a test on `sd_config > 0`. A REML fit at the boundary
returns values like 6e-17 rather than an exact zero, which passes `> 0` while leaving the
deviations with no significant figure to carry: `se` recovers anyway, because it is proportional
either way, while every estimate stays tied at the grand mean. The rule then separates nothing
and the failure is silent. So past `R > 1000` the effects are recomputed from the deviations
instead, which is exact and better conditioned -- at the cost of the `se` being the analytic form
there, since a boundary fit offers nothing to calibrate against.

**3. The estimate has to survive resampling.** This is the one that changed a conclusion. The
reference `sd_config` becomes a **fixed prior with no allowance anywhere for its own
uncertainty**, and it is fitted on 20-odd configurations with one or two rows each -- so it can
be large, confident and meaningless.

*What is "stable" here, and what the gate guards.* The quantity whose stability is in question is
**one number**: the `sd_config` estimated off the reference population. It is not the search's
stability and not a configuration's. It matters because that number is then **pinned** -- handed
to the fit as a known constant -- and the fit has no way to express doubt about a constant. Every
BLUP, every `se`, and therefore the whole in-play set is computed as though it were certain. The
gate guards the pin: drop a third of the reference configurations at random, refit, and see what
comes back; 25 times. If the answer is the boundary -- `sd_config` = 0, "I cannot tell these
configurations apart at all" -- in more than a quarter of the refits, the non-zero value was
carried by a handful of configurations rather than by the population, and it is refused. The run
then keeps its own fitted estimate, which is the pre-2026-11-03 behaviour. **Stability = does
this number survive dropping a third of the evidence behind it.**

*Why the counts below are in the twenties, and whether they stay there.* Those are the **seven
stored runs**, all of which predate `reference_every`. On them the reference population is the
one-time `n_random_init` = 25 cohort, frozen the moment Phase 3 ended -- and `source` was not
recorded before 2026-11-03, so the set is reconstructed from arrival order. Measured sizes, in
configurations that produced a usable score: 20, 20, 24, 24, 25, 25, 25. Going forward the
population **does** grow without bound: `reference_every` = 20 draws a fresh
surrogate-independent configuration one iteration in twenty, forever. But it grows slowly in the
unit that matters -- about one configuration per 20 evaluations, so roughly 500 evaluations to
double a 25-configuration population -- which is why growth and the gate were added together.
Growth is what eventually makes the gate moot; until then the gate is what protects a run whose
reference population is still thin. Neither helps the runs below retroactively.

> **Measured**, refitting each reference population on two thirds of its configurations, 25 draws:
>
> | run | reference `sd_config` | IQR under subsampling | Q3/Q1 | collapses to 0 | usable rows / configs |
> |---|---|---|---|---|---|
> | `34e15bd03838` | 0.0560 | 0.0561–0.0657 | 1.2 | 0/25 | 86 / 25 |
> | `f92e28ec6415` | 0.0251 | 0.0209–0.0299 | 1.4 | 0/25 | 50 / 25 |
> | `9a2664962c01` | 0.0246 | 0.0207–0.0284 | 1.4 | 1/25 | 86 / 24 |
> | `46b0d4b39c79` | 0.0273 | 0.0156–0.0333 | 2.1 | 4/25 | 32 / 20 |
> | `e7517df121d5` | **0.0517** | **0.0265–0.0970** | **3.7** | 5/25 | **37 / 24** |
> | `80cfbd7c27e9` | 0.0130 | 0.0053–0.0188 | 3.5 | 6/25 | 169 / 25 |
>
> `e7517df121d5` has 1.5 usable rows per reference configuration, and pinning its 0.0517 resolves
> its field to **2 or 3 configurations out of 300** after 2964 evaluations.

`reference_stability_draws` (25) runs that resampling at every fit and refuses the estimate when
more than a quarter of the draws collapse to the boundary. A row-count floor would not have done
it: `80cfbd7c27e9` has 6.8 usable rows per configuration and is the least stable of the seven by
that count, because its 25 configurations are spread over 24 trials and few of them share one.
Resampling measures whether the answer depends on which configurations were drawn; a row count
is only a proxy for that.

Two properties the gate needs and did not get for free. The subsample size was
`max(reference_min_configs, 2n/3)`, which with 24 configurations and a floor of 20 can drop only
four, so the probe was weaker than the threshold was calibrated against; it is now
`max(12, ceil(2n/3))`, independent of the acceptance floor. And the draws came from the ambient
RNG, which made the verdict a coin flip **re-rolled at every call** — the same evidence could
pin 0.0517 on one iteration and fall back on the next, the very thrashing §7 is for. The draws
are now seeded from a hash of the reference configurations, and the caller's stream is restored
on exit.

**What it catches, measured.** On this store the gate refuses only `7d75f214584a`, whose
reference fit returns a boundary 0 before the gate is reached at all, plus the first two
checkpoints of `80cfbd7c27e9` while its population is still thin. The other runs sit at 0–24%
collapse, under the quarter, so **the gate is a guard against a reference population that is
genuinely too thin rather than a discriminator among these seven** — and `e7517df121d5`'s 0.0517
is pinned with quartiles spanning a factor of 3.7 that nothing downstream allows for. §15 carries
that as an open item; `LESSONS.md` #44 carries how the earlier, stronger claim about this gate
came to be made.

### The other half: a threshold the search cannot move

A prior estimated off-trajectory stops the *unit* shrinking. It does not supply a notion of how
large a difference is worth chasing. That is `resolution_delta` in §7: a fixed, externally chosen
Δ — the difference in predictive ability actually worth resolving — which does not move when the
population homogenizes. **Measured, the 0.02 this section suggested before it was built is far too
large:** at Δ = 0.02 three fifths of the field is in play. 0.002 is where the trade sits, and §7
has the grid — and the honest reading of that grid, which is that Δ buys set size here and very
little stability. `BACKGROUND.md` ("setting replication from the variance the run itself measures")
records the related design — replication targets driven by Δ rather than by `se` — which is still
not built.

---

## 12. Other search behaviours worth knowing

A comprehensive list of the remaining rules that change what gets evaluated. Each is a fact
about the current code, not a recommendation.

1. **The store is global; a run reads only its own slice.** Universe, scheme and build filters
   compose (§4.1). A configuration with 50 evaluations in another domain is "untried" here.
2. **Build invalidation is per-change, not by version cutoff.** `BUILD_CHANGES`
   (`R/optimizer.R:45`) holds one predicate per result-changing build, so only configurations a
   change actually touched lose their old rows.
3. **The parent pool has no replication requirement.** `get_parents` ranks the in-play set by
   UCB regardless of `n_ok`, so a configuration with one usable evaluation can be a parent.
   Shrinkage is the only guard -- and it is a real one: measured, the median `n_ok` of the UCB
   top 8 is 4.5 to 5 at `z <= 0.75`, because a singleton's BLUP is pulled so far toward the
   grand mean that its inflated `se` does not recover the gap. It falls to 2 only at `z >= 1.5`.
4. **The EI reference `best` has no replication requirement either**, so one lucky
   single-evaluation configuration raises the bar for every candidate.
5. **`incumbent_min_reps` does not affect the search.** `incumbent_config` is called only by
   `write_report` (`R/report.R:119`). The reported incumbent and the configuration the search
   is actually exploiting can differ.
6. **Failures count toward the replication floor** (`base` uses `n_eval`, which includes rows
   with `score = NA`), but not toward `n_ok`. A configuration that fails twice has met the floor
   and, having no finite `mean_score`, can never be a contender either -- so it is never tried
   again (subject only to the §7 exception). The comment on `.replication_backlog` names this
   as deliberate: counting evaluations rather than distinct trials is what makes the tier
   always drain.
7. **The same pair is never evaluated twice** -- the pipeline is deterministic in
   `(config, trial, scheme)`. `duplicate_cells` counts violations and the report prints them.
8. **A configuration that has covered the whole universe leaves the backlog**, and when that is
   true of every contender the report says no further evidence about the leaders is obtainable.
   An unreachable trial in the universe prevents this from ever firing (§6).
9. **Worker offsets, not locks, spread the work.** Seeds, the replication pool and the trial
   backlog are all indexed by `((worker - 1) mod |pool|) + 1`; the contender tier is additionally
   rationed on `(nrow(evals) + worker) %% replicate_every`. Correctness comes from `claim_eval`,
   not from the offsets.
10. **Several knobs in the configuration space are sampled but inert in some branches** --
    `lambda_fixed` is sampled always and used only when `lambda_select = "fixed"`;
    `rkhs_theta` only for `rkhs_gaussian`; `secondary_min` is ignored when
    `primary_only = "yes"`. They still occupy feature columns in the surrogate and still make
    two configurations distinct by hash, so the search can spend evaluations distinguishing
    configurations that compute identically.
11. **The two accession floors are different settings, and only one of them defines the run.**
    `min_trial_acc` (30) decides whether a trial may be *sampled at all* and is one of the ten
    run-defining keys, so changing it starts a new `run_id`. `min_test_acc` (10) decides whether
    a score exists and whether it is usable (§3), and is **not** run-defining -- nor is
    `weight_cap`. Changing either re-analyses the existing rows of the live run in place rather
    than starting a fresh run, which is what makes a re-weighting retroactive but also means the
    store holds rows gathered under earlier values of both. §13 is why the line falls there.
12. **Marker density is not in the search space.** It is fixed by
    `settings$dosage_budget_bytes` at parse time (`LESSONS.md` #16), so rows from runs at
    different budgets are not comparable on that axis; `dosage_budget` is stored per row.
13. **The universe is pinned by trial id, once per run**, so a trial renamed mid-run cannot
    silently leave or join it. That protection depends on the `runs` record surviving, which it
    did not for a while (`LESSONS.md` #32).
14. **`sd_resid` confounds measurement noise with `config x trial` interaction**, because the
    interaction is not identifiable without duplicate cells (§4.3). Any statement of the form
    "the evaluation noise is X" is an upper bound.
15. **The fit's fallback is silent except for `estimator_note`.** A pooled estimator line in
    `report.md` with no note means a choice; with a note it means the fit declined, and then
    contender replication is off entirely and the parent pool falls back to the `n_elites`
    floor (§8.5).
16. **`sd_config_source` changes the `se` estimator, not only the prior.** With a reference
    estimate available the configuration effects are re-shrunk analytically; without one they
    come from lme4's `postVar`. The analytic form runs a few percent low, so a reference run and
    a fitted run are not comparable to that precision (§11).
17. **`var_comps` gains two entries only when a reference estimate was used.** `sd_config` is
    always the fit's own REML estimate so the number never changes meaning; `sd_config_ref` and
    `sd_config_used` appear beside it when the pin was applied. The *absence* of `sd_config_ref`
    in `report.md` is therefore the report of a silent fallback to the fitted value.
18. **The in-play set has a cap for replication and a different cap for parenthood.** The first 8
    by UCB earn extra evaluations; the first `n_parents` = 32 are parents. A configuration ranked
    9th to 32nd breeds but is not deepened.
19. **Under the default `run_identity = "target"`, none of the settings that decide what is in
    play is run-defining.** `contender_rule`, `resolution_delta`, `sd_config_source`,
    `reference_every` and -- under `"target"` -- `contender_z` are all policy, so the selection
    rule is re-tunable on a live run. That is the point of §13, and also the thing to be careful
    about: two stretches of one run's history can have been selected under different rules with
    nothing in the row to say so.
20. **The persistent reference rate spends the acquisition budget, not the replication budget.**
    Phase 3b sits after replication and before acquisition, so `reference_every = 20` costs 5%
    of the steps that would otherwise have gone to EI, and 0% of the `config_replication` floor.
21. **Reference replication is rationed; the `config_replication` floor is not.** A reference
    configuration short of `reference_replication` joins `extra`, the tier `replicate_every`
    rations, so a large `reference_replication` cannot starve the search (§11).
22. **Reference replication requires a recorded `source`; the reference `sd_config` estimate
    does not.** `.reference_hashes(tagged_only = TRUE)` for the first, `FALSE` for the second.
    Spending an evaluation is an action and rests on a fact; estimating a variance from rows
    already in hand is retrospective, and on a store written before 2026-11-03 the arrival-order
    inference is the only thing available. That inference also fails in the direction that
    matters for spending: on a small or untagged store "the first `n_random_init`
    configurations" is most of the store.

---

## 13. What a run is: target or experiment

```r
.RUN_SIGNATURE_KEYS_TARGET <- c("optimize_scheme", "target_domain", "focal_trait_db_id",
                                "min_trial_acc", "simulate")
.RUN_SIGNATURE_KEYS_POLICY <- c("trial_replication", "config_replication",
                                "replicate_every", "contender_z", "n_random_init")

.run_signature_keys <- function(settings)
  if (identical(settings$run_identity %||% "target", "experiment")) .RUN_SIGNATURE_KEYS
  else .RUN_SIGNATURE_KEYS_TARGET

run_id_for <- function(settings, build)
  substr(rlang::hash(c(settings[.run_signature_keys(settings)], list(build = build))), 1, 12)
```
(`R/store.R`.) Same keys, same build, same `run_id`, and the run resumes — inheriting the trial
universe pinned when it was first recorded.

**Apply the test "does this change what a stored row *means*" and the ten keys split in two.**

| | keys | why |
|---|---|---|
| **The target** | `optimize_scheme`, `target_domain`, `focal_trait_db_id`, `min_trial_acc`, `simulate` | Each changes what is being optimized, or what a row *is*. A CV0 score and a CV00 score are different quantities; `min_trial_acc` changes which trials are in the universe at all (§6); `simulate` makes the row synthetic. Mixing rows across any of these is a category error. |
| **Search policy** | `trial_replication`, `config_replication`, `replicate_every`, `contender_z`, `n_random_init` | Each changes the *trajectory* the search takes, not the interpretation of any row it produces. A row written under `config_replication = 2` says exactly what the same row would say under 3: this configuration scored this correlation on this trial under this scheme. |

`settings$run_identity` chooses which half defines a run.

- **`"target"`** (the default) — the five target keys. Policy can be retuned on a live run
  without starting a new one, so evidence accumulates on one target while the method improves.
  This is the reading for *"keep nudging toward a better configuration on this domain."*
- **`"experiment"`** — all ten. Each policy is its own run, so a trajectory is always the
  product of a single policy and a measured churn series is never a blend of two schedules.
  This is the reading for *"compare search policies against each other."*

`run_identity` is in neither set, and cannot be: it selects the signature, so including it would
make a run's identity depend on how that identity is computed.

### What a changed `run_id` does and does not cost

This is worth stating precisely, because the intuition is wrong in a way that matters.

**It does not make the old rows unreachable.** Nothing in selection filters on `run_id`.
`choose_config` composes exactly three filters — `filter_evals_to_universe`,
`filter_evals_to_scheme`, `filter_evals_to_build` (`R/optimizer.R`) — and the store is one
archive shared by every run. A new run over the same domain, scheme and build reads every row
the old run wrote. **No migration, and no new store, is involved in switching `run_identity`.**

What a new `run_id` actually costs is two things:

1. **The trial universe is re-resolved.** `pin_trial_universe` finds no matching row in `runs`,
   so the leader calls `resolve_trial_universe` against the *live* catalogue and pins the result
   (`R/data_access.R`). If the catalogue still resolves the same trial ids, the universe is
   identical and nothing is scoped out. If a trial was renamed, dropped, or has become
   infeasible, the new universe is smaller — and `filter_evals_to_universe` *is* a gate, so rows
   on a trial that left the universe stop counting. This is the real risk, and it is the one
   that bit before (`LESSONS.md` on the replicate livelock: a mid-run rename).
2. **The run record splits.** Any analysis that groups by `run_id` — `dev/sweep_contender_z.R`,
   `dev/sweep_contender_rule.R`, `dev/measure_churn.R` — sees one target's history as two runs.

Switching `run_identity` rotates `run_id` **once**, because the two hashes differ. After that
the chosen reading is stable: under `"target"`, a policy change no longer rotates it at all.

A genuine wrinkle under `"target"`: `record_run` is `INSERT OR IGNORE`, so the `settings_json`
stored against a run is whichever policy recorded it *first*. The report always prints the
**live** policy, so nothing displayed is stale, but the stored snapshot is a record of the run's
first policy rather than its current one. `explain_run_mismatch` names policy differences
explicitly as *"(policy, not identifying)"* so that a resumed run never looks like it matched by
accident.

### Build 0.9.0, and why it retires no rows

`build` is hashed into `run_id` alongside the settings, so a bump rotates the id exactly as a
target-key change does, and the two costs above are its costs. What a bump does *not* do by
itself is retire evidence: that is `BUILD_CHANGES`, one predicate per result-changing build,
and `filter_evals_to_build` drops a row only when some entry's `affects(cfg)` is TRUE and the
row is older than that entry.

**0.9.0 has no entry, deliberately.** The test is whether a change alters what a configuration
*computes on a trial*. Everything the 0.9.0 work touched is downstream of that:

| change | what it alters |
|---|---|
| racing bar, `contender_z`, `resolution_delta` | which configurations are replicated and bred |
| one-scale interval (§7) | how a stored score is *summarized*, not what it is |
| reference-population `sd_config` and its pin (§11) | the prior in an aggregate fit |
| `run_identity` (above) | which rows share a provenance key |
| surrogate case weights | which configurations are proposed |
| `get_parents` cap and floor | which configurations are proposed |

None of those appears anywhere in `run_pipeline()`. A score stored under 0.8.8 is the same
number this build would compute, so every row stays comparable and the whole 49,833-row history
carries forward.

**The one exception, and the judgement call.** `min_test_acc` rose from 5 to 10 in the same
batch, and that *is* a pipeline floor: `.best_panel_index()` admits a panel only at
`f_cov >= min_test`, so raising it can change which panel a single-panel kernel is built from,
and `run_pipeline` now declines a task with 5-9 test accessions outright. An entry predicated
like 0.8.8's — `vanRaden_single`, `rkhs_gaussian`, `best_single_project` — would retire **41,334
of 49,833 rows (83%)**. Against that: the floor changes the panel only when the winning panel
covered 5-9 focal accessions *and* another panel clears both floors at the new threshold; when
none clears them the fallback is `which.max(f_cov)`, which is unaffected and usually returns the
same panel. And the 2,931 stored `ok` rows with `n_test` in 5-9 were already being declined at
aggregation by the same floor, so the search's view of them does not change at all. Retiring
83% of the evidence to repair that was judged the wrong trade. The cost of the judgement is that
a 0.9.0 store may hold a few single-panel rows built on a panel this build would not have
chosen; it is recorded here, and in `BUILD_CHANGES` itself, rather than silently.

### The remaining misfits

`n_random_init` was the clearest misfit in the policy half even before the switch existed,
because the thing it gates is already spent: Phase 3 fires on `n_scored < n_random_init`, a
one-time test against a cumulative count, so raising it on a mature run changes nothing at all.
"Restart exploration on a mature target" needs a persistent exploration *rate*, which is what
`reference_every` now is (§5 Phase 3b, §11).

Four settings that decide what is in play sit outside both signatures: `contender_rule`,
`resolution_delta`, `sd_config_source` and `reference_every`. Under `"target"` so does
`contender_z`, which means the whole selection rule is now re-tunable in place — which is the
point, and also the thing to be careful about, since two stretches of one run's history can have
been selected under different rules with nothing in the row to say so. `evals.source` records
the *phase*, not the policy. If that distinction ever has to be recoverable, the run record is
where it would go.

Two settings sit outside for the same reason read in reverse: `min_test_acc` and `weight_cap`
change how stored rows are *weighted and ranked*, never what any row would have recorded, so
changing them re-analyses the live run in place (§12 item 11). `surrogate_min_obs` is outside it
too, as policy.

---

## 14. What to read in `state/report.md`

The report prints the live value of most of the above (`R/report.R:175-225`):

- `config score estimator: blup (sd_trial …, sd_config …, sd_resid … at unit weight)` -- if it
  says `pooled`, read the note after it; contender replication is idle.
- `contenders: N within racing z = 0.25, delta = 0.002` -- **N is the uncapped in-play count**,
  the resolution series of §9, and the line names the rule and floor it was computed under so a
  report is readable without knowing which settings produced it. `N = 1` means no other
  configuration's claim reaches the bar *under this rule, at this `z`, under this fit*; the
  report calls that "settled", which is shorthand for "nothing else is still in play at this
  confidence level", not a demonstration that the leader is best. `N` equal to the whole scored
  population is the §11 collapse.
- `shrunk on sd_config X from the reference population`, on the estimator line -- which
  `sd_config` the shrinkage used. `shrunk on the FITTED sd_config (no reference estimate
  available)` is the silent fallback of §11 made visible.
- `replication backlog: A base, B contender; C in flight` -- if A is rarely 0, the contender
  tier is being preempted (§5).
- `Contender methods (top 8 by UCB, racing rule at contender_z = 0.25)` -- the capped list, which
  is the replication tier and the first 8 of the parent pool; compare its membership between reports
  to see the turnover series of §9 directly.

The three churn series themselves are not computed by anything in the repo; they are derived
from the store, as in the measurements throughout this document.

---

## 15. What is not known

- **How much the configurations differ from each other, on the 24-trial targets.** `n*` exceeds
  the universe there, so the bound of §4.5 leaves every configuration's posterior close to the
  prior however the budget is spent. `80cfbd7c27e9` nonetheless ended with only 8 configurations
  in play under the uncapped weight, so a top *group* did stand out from the field; whether its
  eight members differ from one another is the part the data cannot answer, and no amount of
  further replication on that universe would settle it.
- **How much of `sd_resid` is `config x trial` interaction.** Not identifiable without
  duplicate cells, which the design counts as a defect.
- **Whether the greedy-EI collapse of §8.3-8.4 is reachable from a healthy run, or only from one
  that starts with a weak signal.** Measured at one point per run (the final state), not traced
  over time.
- **Whether a shared trial panel would fix it**, as `BACKGROUND.md` argues it should. Not built,
  not tested.
- **How much the pinned reference `sd_config` should be discounted for its own uncertainty.**
  It enters the fit as a known constant. The stability gate (§11) only refuses an estimate that
  *collapses* under resampling, and measured on this store that is a weak test: it passes
  `e7517df121d5`, whose draws span 0.0265-0.0970, a factor of 3.7 across the quartiles, against
  1.2-1.4 on the three runs with a comfortable reference population. A run with a wide but
  non-collapsing reference estimate therefore gets a confident prior it has not earned, and the
  effect is largest exactly where it matters most -- the runs whose fitted `sd_config` is 0 and
  which have nothing else to fall back on. Integrating over that uncertainty (a scaled-inverse-
  chi-square prior on the component, or widening the pin by the resampling spread) is the obvious
  next move and is not built; a dispersion criterion in place of the collapse count would be
  tuning a threshold on seven runs, which is why there is none.

---

## See also

- `DESIGN.md` §3 -- the same loop in summary, inside the architecture.
- `BACKGROUND.md` -- why this family of optimizer, and the deferred design for setting
  replication from the measured variances.
- `LESSONS.md` #19 (trial heterogeneity), #30 (the fit that was thrown away), #33 (churn),
  #34 (one floor, not two), #35 (the surrogate's unweighted rows), #36 (the weight cap),
  #37 (the EI reference point, and the exploration term it cost), #38 (selection provenance),
  #39 (the coupled fit threshold).
- `../dev/measure_churn.R` -- every measured number in sections 9-11, re-derived from the
  store through the optimizer's own `aggregate_scores()`, `.contenders()` and candidate
  generators. Run it rather than trusting the figures quoted here: they were taken on the
  store of 2026-10-30 and a later store will differ.
- `../tools/report_replication.R` -- replication as the store currently has it.
