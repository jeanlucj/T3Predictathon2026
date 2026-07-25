# Evaluation guide — stepping through the optimizer

This is the document for **evaluating and validating** the optimizer from the RStudio console -- for me now and in six months. `README.md` is for running it; `DESIGN.md` is the static architecture (every function, where outputs go); `BACKGROUND.md` is why it is built this way and the challenges it must handle. This file is the **runbook**: it lists the modules, gives a fast→slow order to check them in, and walks each one step by step -- one console line at a time -- saying what each step should return, what to eyeball, and the subtle failure signatures to hunt for.

The system is genuinely complex -- the task and the messy T3 data layer make that unavoidable -- so the goal is that **working through this document is faster and more thorough than re-reading all the code**, and that it maximizes the chance of catching the bugs that make the optimizer *look* like it works without delivering on its objective.

The objective, to keep in view the whole time: discover one pipeline configuration whose **mean Pearson correlation between predicted and observed per-accession BLUEs, over randomly sampled T3 grain-yield trials**, beats every one of the five submitted algorithms. A bug that silently hides data or leaks the answer will still produce a plausible-looking number -- that is what we are guarding against.

------------------------------------------------------------------------

## 1. The modules, in one screen

Eleven `R/` modules, one entry point, one settings file. The **group** column is the name you pass to `arm_evaluation()` (§2) to `debug()` that module.

| Module (`R/…`) | Owns | Off/online | Group |
|----|----|----|----|
| `config_space.R` | the six-subtask genome: sample / encode / crossover / mutate | offline | `genome` |
| `seeds.R` | the five submissions as starting configs | offline | `genome` |
| `optimizer.R` | candidate generation + acquisition + phase logic + incumbent | offline | `engine` |
| `surrogate.R` | bagged-rpart RF surrogate: mean, uncertainty, Expected Improvement | offline | `engine` |
| `store.R` | SQLite results store (resumable state) | offline | `store` |
| `report.R` | learning curve, incumbent, method importance, failure log | offline | `store` |
| `evaluate.R` | run a config on a trial → score; the SIMULATE objective | offline\* | `scoring` |
| `conditions.R` | `infeasible` / `fatal` / `sample_failed` signals + funnel | offline | — |
| `data_access.R` | BrAPI + caching: trials, phenotypes, accessions, projects, dosages | **online** | `data` |
| `pipeline.R` | the parameterized six-subtask pipeline + dispatchers + `mask_cv` | **online** | `subtaskA…F`, `flow` |
| `diagnostics.R` | canary oracle (anchor / calibrate / check / coverage) + `diagnose_trial` | **online** | `diagnostics` |

\* `evaluate.R` is offline in `simulate = TRUE`; in real mode it calls `run_pipeline()`, which is online.

The six subtasks inside `pipeline.R` (the recombinable "genome"): **A** select training trials, **B** preprocess phenotypes → per-accession targets, **C** select genotyping data, **D** build relationship/kernel, **E** train model, **F** predict. Each is its own group (`subtaskA` … `subtaskF`) so you can arm exactly one.

------------------------------------------------------------------------

## 2. How to evaluate: the plan

Three principles, in priority order:

1.  **Offline before online.** Prove the math with `simulate = TRUE` (and the test suite, §7) before spending a single download. Levels L1–L4 are offline.
2.  **Fast before slow.** Follow `EVAL_ORDER` (printed by `eval_groups()`); it runs the instant offline modules first and the network/heavy ones last, so a cheap bug is caught before an expensive one. Within the live layer, pick the *smallest* work first (§5).
3.  **One group armed at a time.** `arm_evaluation("subtaskB")`, run the pipeline, step through *only* subtask B, `disarm_evaluation()`. Arming everything at once makes the debugger unusable.

### The three tools

| Tool (in `R/evaluation.R`) | Use it to | When |
|----|----|----|
| `arm_evaluation(group)` / `disarm_evaluation()` | drop into `debug()` for one module and watch the Global Environment change line by line | stepping through logic |
| `peek(x)` | print a one-line health summary of the object flowing between subtasks (shape, NA counts, degeneracy, rowname overlap) and pass it through | inspecting an intermediate |
| `diagnose_trial(id, …)` / `canary_anchor()` / `calibrate_canary_trials()` | compare the pipeline's own counts to an **independent** re-derivation / ground truth | confirming a real result is real |

`eval_groups()` prints the whole menu. The debugger keys, once: **`n`** next line, **`s`** step into a call, **`c`** finish this function (run to its return), **`Q`** quit the debugger. **Always `disarm_evaluation()` before a background run** — otherwise every pipeline call blocks on the debugger.

`peek()` recognizes the five shapes that carry the pipeline's state and flags the signatures that otherwise pass silently: - **named numeric** (`targets`, `obs`, `pred`) → length, distinct count, NA count, range; warns if all finite values are identical. - **matrix** (a dosage matrix, `K`) → dims, square?, symmetric?, diagonal range; pass `accessions=` to get **rowname overlap** (the synonym-mismatch signal). - **list of matrices** (`dosage_list`) → per-project dims and overlap. - **tibble** (`train_obs`) → rows×cols and **NA count per column** (flags `ALL NA`).

------------------------------------------------------------------------

## 3. Bootstrap (paste once per session)

``` r
setwd("scripts/Analysis_Claude/optimizer")   # if not already the working dir
library(tidyverse)
source("settings.R")
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
source("run_optimizer.R")                     # defines run_optimizer()/optimizer_step()
                                              # (its auto-run is guarded, so sourcing is safe)
s <- optimizer_settings()
eval_groups()                                 # confirm the tooling loaded
```

For the **online** levels (L5 on) also open the live connection once and reuse it. The T3
server now requires login, so `t3_connect()` constructs the connection **and logs it in**
from the `T3_USERNAME` / `T3_PASSWORD` environment credentials (copy `.Renviron.example` to
`.Renviron`, fill it in, restart R):

``` r
s    <- modifyList(optimizer_settings(), list(simulate = FALSE))
conn <- t3_connect(s)          # createBrAPIConnection + t3_login from .Renviron creds
```

> **Auth.** If a call comes back `Unauthorized (HTTP 401)` (e.g. the token expired mid-run),
> `.brapi_try()` re-logs in from the same credentials and retries -- no action needed. If the
> re-login itself fails it says so **loudly** (`... re-login FAILED: <reason>`), and if the
> credentials are simply missing it fails **fast** with a message telling you to fill in
> `.Renviron` and **RESTART R** (it is read only at startup, so editing it mid-session does
> nothing -- check with `Sys.getenv("T3_USERNAME")`). Prefer `t3_connect()` over a raw
> `createBrAPIConnection()`: it logs in up front, so missing creds error *immediately* rather
> than after hours of silent 401s.

> **Caveat.** `conn$vcf_archived()` can prompt interactively to pick a file and will hang a non-interactive run. When you script a real-mode probe with `Rscript`, redirect stdin: `Rscript probe.R </dev/null`. Interactively in RStudio you will simply see the prompt.

Durable state lives in two dirs and is regenerable: `state/` (the SQLite store + `report.md`; `touch state/STOP` halts a run) and `cache/` (downloaded data). See §8.

> **Cache backup (remote server).** With `remote_server = TRUE` the cache stays on the work
> disk and is backed up to `$OPTIMIZER_PATH/cache`. `run_optimizer()` reconciles and backs it
> up automatically, but you can also drive it by hand outside the loop -- e.g. before a
> `check_canaries()` / `calibrate_canary_trials()` session that isn't launched through
> `run_optimizer()`:
> ``` r
> restore_cache_from_backup(s)   # warm the WORK cache from the backup (additive: fills gaps,
>                                # never deletes; run it whenever, empty cache or not)
> sync_cache_to_backup(s)        # back the work cache up to durable storage now
> ```
> Both are no-ops when `cache_backup_dir` is unset (local mode) or `rsync` is missing.

------------------------------------------------------------------------

## 4. Level-by-level walkthrough

Each level: **arm** the group, run the numbered console lines (each assigns a Global-Environment variable), and check the three annotations — **→ returns** (shape you should see), **🔍 eyeball** (the `peek()`/inspection and what healthy looks like), **🚩 red flag** (the subtle failure to hunt), plus the italic *Assumption:* the code is making at that step (test **it**, not just the output). `disarm_evaluation()` when you finish a level.

### L1 — `genome` (offline, instant)

``` r
arm_evaluation("genome")
cfg  <- sample_config()                        # step with n/c to watch each block fill
cfgs <- seed_configs()                         # the 5 submissions as configs
kid  <- crossover(cfgs$Prediction1, cfgs$Prediction5)
mut  <- mutate_config(cfg, n_subtasks = 1)
disarm_evaluation()
```

- **→ returns** `cfg` is a named list with one `*.method` per subtask plus its parameters; `cfgs` has 5 named entries; `kid`/`mut` are configs.
- **🔍 eyeball** `config_hash(cfg)` is stable across repeated calls; `str(cfg)` shows every subtask present; each block in `kid` matches *one* parent (`kid$kernel.method %in% c(cfgs$Prediction1$kernel.method, cfgs$Prediction5$kernel.method)`).
- **🚩 red flag** a parameter that does not apply to its method is non-`NA` (should be `NA`), or a required applicable parameter is `NA`; `mut` identical to `cfg` (mutation must change ≥1 block).
- *Assumption:* `config_space.R::SUBTASKS` is the single source of truth — every method named here must have a dispatcher branch in `pipeline.R` (§10).

### L2 — `engine` (offline)

Build a tiny store of synthetic evaluations, then watch the decision engine. `choose_config()` runs in phases, so seed the store the way the loop itself does — otherwise it just hands back the five submissions forever (Phase 1).

``` r
s   <- modifyList(optimizer_settings(), list(simulate = TRUE))
con <- open_store(tempfile(fileext = ".sqlite"))

# On a FRESH store, choose_config()'s Phase 1 returns the five submissions one by
# one and keeps doing so until EACH seed has a row in the store -- so step it once
# now to see `source = "seed:Prediction1"`, then seed them for real:
choose_config(con, s)$source                   # -> "seed:Prediction1" (nothing stored yet)
seeds <- seed_configs()
for (r in 1:3) {                               # a few trials each -> reps for the incumbent
  tr <- sample_trial(s)
  for (nm in names(seeds)) for (sc in s$schemes) {
    ev <- evaluate_config_on_trial(seeds[[nm]], tr, sc, s)
    store_eval(con, seeds[[nm]], tr$id, sc, ev$score, ev$n_test,
               ev$status, ev$reason, ev$detail %||% NA, ev$seconds)
  }
}
choose_config(con, s)$source                   # -> "random_init" (Phase 1 cleared, <n_random_init scored)

# Now fill with random configs until >= n_random_init are scored, so the surrogate
# can take over (Phase 3):
for (i in 1:30) { tr <- sample_trial(s); cfg <- sample_config()
  for (sc in s$schemes) { ev <- evaluate_config_on_trial(cfg, tr, sc, s)
    store_eval(con, cfg, tr$id, sc, ev$score, ev$n_test, ev$status, ev$reason, ev$detail %||% NA, ev$seconds) } }

arm_evaluation("engine")
ch <- choose_config(con, s)                    # now steps INTO the acquisition phase
disarm_evaluation()
```

- **→ returns** `ch` is `list(cfg, source, ei)`; `source` is one of `seed` / `random_init` / `reeval` / `acquisition`.
- **🔍 eyeball** with few evals `source` is `seed` then `random_init`; only past `n_random_init` does it become `acquisition`. Inside acquisition, `expected_improvement()` is always ≥ 0 and rises with predicted mean.
- **🚩 red flag** the surrogate predicts a flat score for every candidate (no signal learned); acquisition picks the same config forever (EI collapsed); `incumbent_config()` returns a config with fewer than `incumbent_min_reps` reps.
- *Assumption:* the surrogate trains on the **per-config mean** score (`aggregate_scores`), i.e. the best config for a *random* trial is the one with the best average — not tuned to any single trial.

### L3 — `store` (offline)

``` r
arm_evaluation("store")
ev  <- read_evals(con)                          # the whole table as a tibble
write_report(con, s)                            # writes state/report.md
disarm_evaluation()
```

- **→ returns** `ev` has columns `config_hash, config_json, trial_id, scheme, score, n_test, status, reason, detail, seconds, ts`; failures are rows too (with `NA` score).
- **🔍 eyeball** round-trip a config: `config_from_json(config_to_json(cfg))` equals `cfg` **including `NA` params**; `report.md` renders incumbent + method importance + failure log.
- **🚩 red flag** `NA` scores dropped instead of stored (the failure log is how we see data-hiding bugs); JSON round-trip changes a parameter's type.
- *Assumption:* the entire optimizer state is reconstructable from this one table, so a killed run resumes with no loss.

### L4 — `scoring` (offline; the CV0/CV00 crux)

``` r
arm_evaluation("scoring")
x  <- setNames(rnorm(20), paste0("g", 1:20))
score_predictions(x, 2 * x + 1)                 # must be +1 (cor of x with a*x+b, a>0)
tr <- sample_trial(s); cfg <- sample_config()
evaluate_config_on_trial(cfg, tr, "CV00", s)    # step through the sim objective
disarm_evaluation()
```

- **→ returns** `score_predictions` returns a scalar in [-1, 1] or `NA` (with a reason) when \<5 accessions overlap or predictions are constant.
- **🔍 eyeball** `cor(x, x) == 1`; only shared accession *names* join; `mask_cv` under `"CV00"` drops **exactly** the focal-trial accessions that `"CV0"` keeps (`setdiff(names(cv0_targets), names(cv00_targets))` = the focal rows).
- **🚩 red flag** a suspiciously **high** score (leakage — the focal answer reached the training side); **CV0 and CV00 scoring identically** (the CV00 mask is a no-op — the whole challenge collapses).
- *Assumption:* this is the exact Predictathon metric (`T3_predictathon_scripts/analysis/30_All_accuracies…Rmd`).

> **Everything below is online.** Switch to the real connection (§3) and follow §5: smallest / cheapest work first.

### L5 — `data` (live BrAPI; start on the CACHED trial)

Trial `10676` is already in `cache/`, so this touches **no** network on the happy path — the cheapest possible first live check.

``` r
arm_evaluation("data")
trial <- build_trial_descriptor("10676", conn, s)   # step: watch the descriptor fill
obs   <- get_observations("10676", conn, s)
acc   <- get_trial_accessions("10676", conn, s)
disarm_evaluation()
peek(obs); peek(trial$accessions)
```

- **→ returns** `trial` is a list with `id, accessions, n_acc, program, location, year, lat/long/elev`; `obs` is a tidy tibble (`study_id, germplasm_name, value, unit_id, rep, block, col, row`).
- **🔍 eyeball** `peek(obs)` — `value` finite, and **`rep`/`block` NOT all `NA`** (they can be legitimately absent, but all-`NA` means BLUE silently degrades to germplasm means downstream — the canonical subtle bug); `trial$lat` finite if the location has coordinates.
- **🚩 red flag** `n_acc` far below the T3 web UI's count for the trial; a camelCase-vs-snake_case column surprise (raw `conn$search` responses are camelCase, `T3BrapiHelpers` outputs are `janitor::clean_names()`).
- *Assumption:* the `obs_<sid>` cache is trait-independent — the focal trait is filtered at *read* time, so changing `focal_trait` reuses it.

### L6 — `subtaskA` select training trials

Subtask A has **three** methods and they share almost no code, so each gets its own level: **L6a** `accession_overlap`, **L6b** `top_k_similar`, **L6c** `same_program`. All three run on the cached trial `10676`, cheapest first. Arm once for all of them:

``` r
arm_evaluation("subtaskA")     # ... disarm_evaluation() after L6c
```

`.find_related()` derives the neighbour list from the **`acc_<sid>` cache**: one wizard query returns the candidate trials, then the germplasm-overlap counts are set intersections of accession lists. So the first call on a *cold* trial fetches `acc_` for every candidate (watch the "Germplasm overlap: candidate trials" progress bar), and every later call — on any trial, at any threshold — reuses those files. The counts are also memoized in RAM for the session, so repeats are instant. If a *repeat* is slow, the `acc_` cache is not being hit.

#### L6a — `accession_overlap` (start with FEW trials)

``` r
ao <- function(pmin, smin, only) modifyList(seed_configs()$Prediction1,
        list(train_select.method = "accession_overlap", train_select.primary_min = pmin,
             train_select.secondary_min = smin, train_select.primary_only = only))
# HIGH overlap threshold first -> few neighbours -> few downloads later:
cfgA   <- ao(20, 12, "yes")                                        # carried into L7-L10
tr_ids <- prim <- select_training_trials(cfgA, trial, conn, s)     # primary tier only
both <- select_training_trials(ao(20, 12, "no"),  trial, conn, s)  # + secondary tier
lo   <- select_training_trials(ao(20,  6, "no"),  trial, conn, s)  # looser secondary
```

- **→ returns** a character vector of study ids (the focal id is removed downstream in `run_pipeline`).
- **🔍 eyeball** a *high* `primary_min` returns a short list; lowering it lengthens it monotonically. The secondary tier is **additive**: `all(prim %in% both)` is `TRUE`, and loosening `secondary_min` only grows it (`all(both %in% lo)`).
- **🚩 red flag** an empty list where the trial obviously has neighbours (a `.find_related` / id-space bug); a list that ignores the threshold; `setequal(prim, both)` — the secondary tier added *nothing*, i.e. the expansion is silently a no-op.
- *Assumption:* the secondary tier pools the germplasm across **all** primary trials and makes **one** `.find_related` call against that combined accession set (via `.find_related(pooled_acc, …, input = "accessions")`). A candidate is judged on its overlap with the pool *as a whole*, so a trial sharing a few accessions with each of several primary trials can qualify even though it clears the threshold against no single primary trial — this is deliberately more inclusive than the old per-primary-trial hop (which itself replaced a hidden `head(primary, 10)` cap). `primary_min` sets how large the primary pool is; `secondary_min` is the overlap threshold against the pooled germplasm.
- *Assumption:* `.find_related` counts overlap in germplasm **name** space (from `acc_<sid>`, the same space every downstream join uses) and returns only trials in `trial_catalog()` — i.e. **only trials that measured the focal trait**. A neighbour that shares germplasm but never measured the trait is correctly absent: it would have counted toward `min_train_trials` while contributing no observations.

#### L6b — `top_k_similar` (the ranking oracle)

``` r
tk <- function(k, sim) modifyList(seed_configs()$Prediction3,
        list(train_select.method = "top_k_similar", train_select.k = k,
             train_select.similarity = sim))
k8  <- select_training_trials(tk(8,  "genomic"), trial, conn, s)
k16 <- select_training_trials(tk(16, "genomic"), trial, conn, s)

# INDEPENDENT re-derivation of the genomic ranking (does not call .trial_similarity):
pool <- .find_related(trial$id, conn, s, 1)                       # the candidate pool
ov   <- purrr::map_int(pool, ~ length(intersect(
          get_trial_accessions(.x, conn, s), trial$accessions))) |> setNames(pool)
identical(sort(k8), sort(names(head(sort(ov, decreasing = TRUE), 8))))   # -> TRUE
```

- **→ returns** at most `k` study ids — `k` is a **cap**, not a guarantee: the pool can be smaller than `k`.
- **🔍 eyeball** the independent check above: under `similarity = "genomic"` the selected trials are exactly the `k` with the largest shared-accession count. Raising `k` is **nested** (`all(k8 %in% k16)`).
- **🚩 red flag** `similarity = "environmental"` when `trial$lat` is **not finite** — the whole environmental branch is skipped, every score stays `0`, and `head()` silently returns an **arbitrary** `k` trials that look like a real answer. Check `is.finite(trial$lat)` *before* trusting an environmental or `both` selection. Also: `k` ids returned but the ranking is uncorrelated with `ov` → the similarity is being computed against the wrong accession set.
- *Assumption:* the candidate pool is `.find_related(focal, …, min_common = 1)` — i.e. **every** focal-trait trial sharing ≥1 germplasm, often hundreds. `.trial_similarity` then calls `get_trial_accessions()` once per candidate — the same `acc_<sid>` files `.find_related` just populated, so the "Similarity: candidate trials" bar should fly. If it crawls, the `acc_` cache is being missed.

#### L6c — `same_program` (the fallback trap)

``` r
sp <- function(cap) modifyList(seed_configs()$Prediction5,
        list(train_select.method = "same_program", train_select.prog_cap = cap))
p40 <- select_training_trials(sp(40), trial, conn, s)
p10 <- select_training_trials(sp(10), trial, conn, s)

# INDEPENDENT check against the catalogue:
cat  <- trial_catalog(conn, s)
mine <- cat |> dplyr::filter(program_name == trial$program, study_db_id != trial$id)
all(p40 %in% as.character(mine$study_db_id)); length(p40) <= 40   # -> TRUE, TRUE
```

- **→ returns** at most `prog_cap` ids, all from the focal trial's breeding program, focal excluded.
- **🔍 eyeball** every returned id is in `mine` (same program) and the cap is respected. Note the program's total size, `nrow(mine)`: if it is ≤ `prog_cap`, the cap branch never fires and `p40` is simply the whole program.
- **🚩 red flag** the `> prog_cap` branch re-filters to trials at the **same location or in the same year**, and if nothing matches it returns **empty** — a big program can collapse to zero training trials, which then reads downstream as an ordinary `too_few_train_trials` infeasibility rather than as this fallback misfiring. If `p10` is empty while `nrow(mine)` is large, that is the trap, not a genuinely unusable trial. Confirm with `dplyr::filter(mine, location_name == trial$location | year == trial$year)`.
- *Assumption:* `same_program` selects on **metadata only** — nothing guarantees the chosen trials share any germplasm with the focal trial, so this method can hand a perfectly large training set with no genomic connection to subtask C/D. A high trial count here is not evidence of a usable training set.

``` r
disarm_evaluation()
```

### L7 — `subtaskB` preprocess phenotypes → targets

``` r
arm_evaluation("subtaskB")
train_obs <- mask_cv(get_observations(tr_ids, conn, s), trial$accessions, "CV0")
targets   <- build_targets(cfgA, train_obs, trial, conn, s)
disarm_evaluation()
peek(train_obs); peek(targets)
```

- **→ returns** `targets` a named numeric, one value per training accession.
- **🔍 eyeball** `peek(train_obs)` — no column is `ALL NA`; `peek(targets)` — `distinct` ≫ 1 and a plausible yield range.
- **🚩 red flag** **targets all identical / degenerate** (`build_targets` collapsed — e.g. every trial centered to the same mean); `rep`/`block` all `NA` so `.blue_per_trial` fell back to plain means for every trial without your noticing.
- *Assumption:* per-trial BLUE `y ~ germ (+rep +block)` only adds `rep`/`block` to the model **when they vary** — so their absence is tolerated, not detected.

### L8 — `subtaskC` select genotyping data (the synonym check)

``` r
arm_evaluation("subtaskC")
dl <- choose_geno_sources(cfgA, names(targets), trial$accessions, conn, s)
disarm_evaluation()
peek(dl, accessions = c(names(targets), trial$accessions))
```

- **→ returns** `dl` a named list with **one matrix per protocol group** (name = the merged project ids, e.g. `"2762+9441"`), each holding that group's **full genotyped population** — every sample in the constituent projects, not just the accessions this trial needs. `nrow()` in the hundreds-to-thousands is correct and expected. An `n_projects` attribute carries how many covering projects we started from.
- **🔍 eyeball** `peek(dl, accessions=…)` — each matrix's **rowname overlap** with the accession set is \> 0 (it will be far smaller than `nrow`, which is the point). `names(dl)` shows which projects got merged: same-panel projects (a protocol and its `V2`/`v2.1` re-registration against another reference genome) should land together.
- **🚩 red flag** a **non-empty dosage matrix with overlap = 0** — the decisive synonym / name-mismatch signature (VCF samples under old preliminary names); `n_projects > 0` but `length(dl) == 0` (found projects, extracted nothing → a download/parse bug, flagged `suspect` by `run_pipeline`). Also: two obviously different panels merged into one group (check `settings$merge_containment`), or a group whose matrix has far fewer accessions than its projects did (over-eager `.prune_redundant`). Cross-check with `diagnose_trial("10676", s, conn)`.
- *Assumption:* a genotyping **project** id (not protocol id) is what `vcf_archived()` accepts. Protocol *ids* are **not** used to group projects — a `V2`/`v2.1` protocol is the same protocol scored against a different reference genome and carries a different id, so grouping is done on **marker overlap** (`settings$merge_containment`, 0.95).

### L9 — `subtaskD` build the kernel

``` r
arm_evaluation("subtaskD")
need <- union(names(targets), trial$accessions)
K <- build_kernel(cfgA, dl, need)
disarm_evaluation()
peek(K)

# Marker QC no longer depends on WHO is needed. Confirm it:
Q <- .qc_markers(dl[[1]], cfgA)                    # QC on the panel's full population
ncol(Q)                                            # same number for any `need`
ncol(.qc_markers(dl[[1]][sample(rownames(dl[[1]]), 10), ], cfgA))   # the OLD way: fewer
```

- **→ returns** `K` a square relationship matrix over the **needed** accessions only (training ∪ focal, intersected with the genotyped), even though QC and the allele frequencies behind it came from the panel's full population.

- **🔍 eyeball** `peek(K)` — **square, symmetric**; clones (identical dosage rows) give identical `K` rows. The two `.qc_markers` calls above should differ substantially: on a 3,707-sample panel, QC from 10 needed accessions keeps \~6,900 markers against the population's \~9,000 — a quarter of the marker set decided by *who was needed*. That disagreement is the bug this level exists to keep fixed.

- **🚩 red flag** `symmetric = FALSE`; `em_combine` returning a matrix one row/col too big (the phantom-variance row must be dropped); fewer than 50 markers surviving QC (`.qc_markers` should have raised `too_few_markers`); a **`NaN` GRM**.

- *Diagonal expectation — read the method first.* The diagonal only diagnoses anything under `vanRaden_single`:

  | method | expected diagonal | why |
  |----|----|----|
  | `vanRaden_single` | **\~1.5–2.0** on real T3 data | VanRaden's diagonal is `1 + F` and wheat lines are highly inbred. **\~1.0 here is suspicious, not reassuring** — the old miscoded kernel produced exactly that (mean 0.95) because `A.mat`'s scaling denominator was wrong. |
  | `em_combine` | **\~1.0, by construction** | `build_kernel` standardises each partial covariance (`g / mean(diag(g))`) before `covariance_combiner()`. A diagonal near 1 here means nothing. |
  | `rkhs_gaussian` | **exactly 1.0** | `exp(-theta·0)` on the diagonal, by definition. |

  So do not read a \~1.0 diagonal as the coding bug unless the config is `vanRaden_single`.

- *Assumption:* `ridge` is added to the diagonal *after* the base kernel; markers are the intersection across projects within a group (largest-project fallback below 50 shared).

- *Assumption (em_combine):* the per-panel GRMs are built over `need ∪ bridge`, where `bridge = .bridge_accessions(dosage_list)` are accessions genotyped in **\>1 protocol group**. Those shared rows are what `covariance_combiner()` stitches the separate GRMs through; the combined GRM is then subset back to `need` (bridges are scaffolding, dropped from the result). Confirm by stepping: `bridge` should be non-empty when panels overlap (real cached panels share 285 such accessions), and `rownames(K)` after the branch is a subset of `need` with no bridge-only accession left. When `bridge` is empty (disjoint panels) this is a no-op and the cross-panel blocks are legitimately 0 — see §6.

- *Assumption:* **dosages are coded `{0,1,2}`.** That assumption did not disappear when `rrBLUP::A.mat()` was dropped — it moved into `.vanraden()`'s `p <- colMeans(X)/2`, and it is the thing to protect. `.vanraden()` now `fatal()`s on a negative entry (`bug_dosage_coding`), so re-encoding `.vcf_to_dosage` to `{-1,0,1}` fails loudly rather than silently producing a meaningless `p`. The original bug was the mirror image: `A.mat` documents `{-1,0,1}`, so fed `{0,1,2}` it computed `freq = p + 0.5`, its internal `MAF = pmin(freq, 1-freq)` went negative for every alt-major marker, and its `min.MAF` filter dropped them (on an all-alt-major panel: *every* marker, giving an all-`NaN` GRM). Do not "simplify" back to `A.mat` without recoding to `X - 1`.

- *Assumption:* because `K_ij` depends on other accessions only through the allele frequency `p`, building it over `need` with population `p` gives **exactly** the `[need, need]` submatrix of the whole population's GRM — the reordering is for cost, not an approximation. `tests/test_subtasks.R` asserts that identity.

### L10 — `subtaskE` + `subtaskF` train and predict

``` r
arm_evaluation(c("subtaskE", "subtaskF"))
geno <- rownames(K); tin <- intersect(names(targets), geno); ten <- intersect(trial$accessions, geno)
fit  <- train_model(cfgA, targets[tin], K, tin, ten, trial, train_obs)
pred <- predict_test(cfgA, fit, K, tin, ten, targets, "CV0", s)
disarm_evaluation()
peek(pred)
```

- **→ returns** `fit` a model object (`kind = "gblup"`, a `u` random effect, a `mu`); `pred` a named numeric over the test accessions.
- **🔍 eyeball** `peek(pred)` — **not constant**, plausible range; under the fallback (too little overlap) it is deliberately `mean(targets)` for everyone.
- **🚩 red flag** predictions constant when overlap was adequate (the model collapsed, not the fallback); blending changing predictions under **CV00** (it must be a no-op there — those phenotypes were masked).
- *Assumption:* `predict_test`'s `cond_expectation` computes `mu + K21 K11⁻¹ u`; the observed-BLUE blend applies **only** under CV0.

### L11 — `flow` end to end (watch the funnel)

``` r
arm_evaluation("flow")               # run_pipeline + optimizer_step only (no double-break)
po <- run_pipeline(cfgA, trial, "CV0", s, conn)   # c through it, watch stage counts
disarm_evaluation()
score_predictions(po$pred, po$obs)
```

- **→ returns** `po = list(pred, obs)`; the score is a plausible correlation (roughly in the range the Predictathon reported, not \~1 and not `NA`).
- **🔍 eyeball** the stage counts stay sizeable end to end — no **cliff** (many accessions in, \~zero out) between the phenotyped and genotyped sides.
- **🚩 red flag** a cliff `run_pipeline` tags `suspect`; a score too high (leakage) or always `NA` (silent infeasibility).
- *Assumption:* an unevaluable (trial, config) raises `infeasible` and is *recorded* (not a crash); only a wiring/settings error raises `fatal` and halts.

### L12 — `diagnostics` (the independent oracle)

``` r
arm_evaluation("diagnostics")
anchor <- canary_anchor(s)                       # counts from the 5 teams' submission CSVs (no optimizer code)
cal    <- calibrate_canary_trials(s, conn, anchor)   # our counts vs the anchor, per trial (LIGHT: no model fits)
print_calibration(cal)
check_canaries(s, conn)                          # coverage oracle: one frozen config per trial
sweep_rich_trials(s, conn)                       # code oracle: every method/branch on 2 rich trials
disarm_evaluation()
```

There are **three** distinct checks here; run the light one first, then pick the oracle that answers your question:

| Tool | What it asks | Cost | Read a failure as |
|---|---|---|---|
| `calibrate_canary_trials` + `canary_anchor` | Does our per-trial **data count** match independent ground truth (the 5 teams' submissions)? | Light — no VCF download (`deep=FALSE`), no model fits | a `divergent` row = a **data-plumbing bug** localized to that trial |
| `check_canaries` | Can the pipeline predict nine known trials under **one frozen config each**? | Heavy — downloads + fits for 9 trials | ambiguous: `error`=bug, but `infeasible`/`constant`=this config × that trial's data (or a stale freeze) — see "Reading a `check_canaries` result" in §9 |
| `sweep_rich_trials` | Does **every method/branch** produce a prediction on **data-rich** trials (10675, 10677)? | Medium — 2 trials' genotypes (cached after) × a fit per variant | on rich data feasibility is a code property, so a `SUSPECT`/`BUG` variant is a **code bug**, not a data verdict |

**Which to run.** Always do the calibrate step — it's cheap and rests on independent ground truth, so it's the trustworthy plumbing check. Then: reach for **`sweep_rich_trials`** when you want "is the *code* correct?" (after a change to the pipeline), because it isolates a broken branch without the data-adequacy ambiguity. Reach for **`check_canaries`** when you want "does this *frozen production config* still work on the trials I care about?" — but only trust its `infeasible`/`constant` verdicts once its configs have been validated by calibration (they have not been fully, so today read those as "data/method," not "bug").

- **→ returns** `cal` a per-trial table with a `divergent` flag; `check_canaries` prints `CANARY ALARM` (hard, strong trial) or a soft warning (weak trial); `sweep_rich_trials` prints `SWEEP CLEAN` or a `SWEEP ALARM` naming suspect variants.
- **🔍 eyeball** every non-`divergent` calibration row means our pipeline's count **agrees** with the independent anchor — positive evidence the plumbing is right; a `SWEEP CLEAN` means every code path predicts on rich data.
- **🚩 red flag** any `divergent` row (our count outside 0.5×–2× the anchor) — a bug localized to that trial; a `SWEEP ALARM` variant — a broken code branch (cross-check with `diagnose_trial()`); a `check_canaries` `error` — a crash bug (but weigh `infeasible`/`constant` against §9's bug-vs-data table).
- *Assumption:* the anchor shares **none** of the optimizer's code, so a bug rarely corrupts it and our derivation identically — agreement is real evidence. `check_canaries`'s frozen `canary_configs()` are a *coverage* set, not a validated-feasible one (freeze only once calibration agrees, human-reviewed); `sweep_rich_trials` sidesteps that by testing on trials rich enough that a failure is the code's fault.

------------------------------------------------------------------------

## 5. Fast → slow within the live layer

Live work is dominated by VCF downloads and mixed-model fits. Order it so the cheap paths clear first:

- **Cached trial first.** L5–L11 above use `10676` (already in `cache/`) so the first full end-to-end pass downloads nothing. Only then move to a fresh trial.
- **Few training trials before many.** In subtask A, start with a *high* `accession_overlap` `primary_min` (fewest neighbours → fewest downloads); lower it only once the high-threshold path is clean. The prompt's rule of thumb — overlap 20 before overlap 3. Turning `primary_only` off adds a secondary lookup over the pooled primary germplasm, so a low `primary_min` (large pool) plus a low `secondary_min` (loose threshold) is the most expensive corner of the whole search space; go there last.
- **Neighbour lookups warm the `acc_` cache; dosages are what actually cost.** `.find_related()` fetches `acc_<sid>` for each candidate trial on a cold cache, then reuses them for every later trial and for `.trial_similarity`. What the bigger training sets still cost you is downstream: more training trials → more genotyping projects → more VCFs in subtask C.
- **Smallest canary before largest.** `calibrate_canary_trials()` / `canary_anchor()` print per-trial accession counts; check the smallest trial before the largest so a bug surfaces on cheap data.
- **Thin markers on the first touch of a big project.** Pass `marker_thin > 1` (e.g. via a config's `geno_select.marker_thin`) so the first parse of a large VCF is cheap; repeat un-thinned once the path is trusted.
- **Watch the progress bars.** The slow loops report progress: `"Load project dosages"` (per project in `choose_geno_sources`), `"Similarity: candidate trials"`, `"Observations from study ids"`, `"Genotyping projects for accessions"`, `"Canary trials"`, `"Calibrate canaries"`. A bar that stalls tells you *which* item is slow.

------------------------------------------------------------------------

## 6. Subtle-bug catalogue

The failure modes that make a broken optimizer look like a working one, each with the check that exposes it. This is the list to internalize — the whole point is noticing something "a little funny" that the code's own checks did not.

| Signature | What it means | Expose it with |
|----|----|----|
| `rep`/`block` (or any column) **all `NA`** | BLUE silently degraded to germplasm means | `peek(train_obs)` NA-per-column |
| `top_k_similar` returns `k` trials on a trial with **no coordinates** | `similarity = "environmental"`/`"both"` skips its branch when `trial$lat` is not finite → every score 0 → an **arbitrary** `k` | L6b: `is.finite(trial$lat)`; compare the selection to the genomic overlap ranking |
| `same_program` returns **empty** on a large program | the `> prog_cap` branch re-filters to same location-or-year and keeps nothing; reads downstream as ordinary `too_few_train_trials` | L6c: compare with `nrow(mine)` (program size) before concluding infeasible |
| **targets all identical / degenerate** | `build_targets` collapsed | `peek(targets)` distinct count |
| **dosage rowname overlap = 0** on a non-empty matrix | synonym / name mismatch | `peek(dl, accessions=…)`; `diagnose_trial` |
| **K not symmetric / clones differ** | kernel construction bug | `peek(K)` |
| **K diagonal \~1.0 under `vanRaden_single`** on inbred wheat | a dosage-coding bug: alt-major markers dropped and/or the scaling denominator wrong. Correct VanRaden here is `1 + F` ≈ 1.5–2. (Only diagnostic for `vanRaden_single` — `em_combine` standardises its diagonal to 1 and RKHS's is 1 by definition.) | L9 diagonal table; `tests/test_subtasks.R` alt-major oracle |
| **marker count after QC changes with the size of the needed set** | QC is being estimated from the needed accessions, not the panel population (10 needed → \~6.9k markers; population → \~9.0k) | the two `.qc_markers` calls in L9 |
| **em_combine cross-panel blocks all \~0** when panels *do* share accessions | bridge accessions (genotyped in \>1 protocol) are being dropped before the combine, so `covariance_combiner` has nothing to stitch through — the value of `em_combine` collapses to block-diagonal | step `.bridge_accessions(dosage_list)` (should be non-empty); a cross-panel block should be non-zero when a bridge exists, exactly 0 when it doesn't |
| **predictions constant → `NA`**, or **score ≈ 1** | model collapse, or leakage | `peek(pred)`; sanity-check the score |
| **CV0 and CV00 scores identical** | `mask_cv` no-op — the CV00 mask isn't excluding focal accessions | the L4 differential check |
| **funnel cliff** (many in, \~0 out), tagged `suspect` | data-hiding bug, not genuine infeasibility | the `flow` funnel; `report.md` ⚠ section |
| **incumbent never beats the seeds** after many iters | search-space / surrogate bug | `report.md` learning curve |

The meta-rule behind all of them: **trust agreement between two independent derivations; distrust a single number.** Our per-trial count vs `canary_anchor`; the pipeline funnel vs `diagnose_trial`'s raw re-derivation. One bug rarely corrupts both the same way — so a match is evidence, and a mismatch localizes.

------------------------------------------------------------------------

## 7. The test suite (offline, deterministic)

### Principle: oracle tests, not snapshot tests

A **snapshot test** asserts "the function returns whatever it currently returns." It re-encodes the implementation, passes even when the code is wrong, and takes as long to audit as the code. We avoid those. Every test has an **independent oracle** — the expected answer is knowable without reading the implementation. Four kinds:

- **Closed-form** — correlation of a vector with itself is 1; a GRM among genetically identical lines has identical rows; centering by trial mean removes a constant added to one trial.
- **Property / invariant** — Expected Improvement is ≥ 0; a relationship matrix is symmetric; two-stage BLUP shrinks effects toward the mean; QC drops a monomorphic marker.
- **Metamorphic / differential** — double the phenotypes → BLUEs double; permute accessions → predictions permute identically; **CV00 excludes exactly the focal-trial phenotypes CV0 keeps**; raising `blend_obs_w` changes CV0 predictions but is a no-op under CV00.
- **Tiny hand-computable** — a 3-marker × 2-sample VCF whose dosage matrix you can write down by hand.

Because the tests *run*, a test that encodes a wrong expectation usually **fails against correct code**, surfacing the disagreement to adjudicate rather than silently lying. Each test file comments its oracle inline.

Unit tests cover the **deterministic, offline** logic — most of the optimizer. They do **not** verify that the live BrAPI responses have the shape `R/data_access.R` assumes; that is what the canary oracle (§9) is for. *Unit tests prove the math; the canary proves the plumbing.*

### Running, and expected outcomes

``` bash
cd scripts/Analysis_Claude/optimizer
Rscript tests/run_all.R          # fast files; one PASS/FAIL summary
Rscript tests/run_all.R --all    # also the slow end-to-end sim loop
Rscript tests/test_subtasks.R    # a single file
```

Each file is self-contained: it sources the subsystem, runs hand-rolled `check()` assertions (no `testthat`, to match project style), prints a pass/fail count, and exits non-zero on failure. `run_all.R` runs each in a fresh process and aggregates.

| Command | Expected |
|----|----|
| `tests/test_config_space.R` | \~5 genome invariants across \~400 sampled/recombined configs → `config_space tests: 8007 passed, 0 failed` (8007 = individual assertions) |
| `tests/test_subtasks.R` | `Tier 1 subtask tests: 144 passed, 0 failed` |
| `tests/run_all.R` | `2/2 test files passed` |
| `tests/test_sim_loop.R` (or `run_all.R --all`) | `PASS: optimizer beats submissions and improves over random search`, exit 0 |

If a count drifts after a change, that is the regression signal — reconcile it before moving on.

### Tier 1 — deterministic core *(implemented: `tests/test_subtasks.R`)*

| Component (function) | Oracle |
|----|----|
| `score_predictions` (`evaluate.R`) | cor(x, a·x+b)=1 for a\>0, =−1 for a\<0; only shared accession names join; \<5 overlap → NA `too_few_overlap`; constant → NA `constant`. |
| `.vcf_to_dosage` (`data_access.R`) — guards the transposition bug | hand 3-marker × 2-sample VCF → **accessions × markers**, `0/0→0,0/1→1,1/1→2`, correct names; sample subset → one row; `thin=2` keeps every 2nd marker. |
| `.qc_markers` (`pipeline.R`) | monomorphic dropped (MAF=0); high-missing dropped; remaining NA imputed (mean vs mean_round); errors below 50 markers. |
| `.merge_markers` (`pipeline.R`) | keeps marker intersection; de-dups repeated accessions; falls back to largest project when \<50 shared; single project as-is. |
| `build_kernel` (`pipeline.R`) | VanRaden K symmetric, clones give identical rows; `ridge` raises the diagonal by exactly `ridge`; RKHS K symmetric, diag 1, off-diag (0,1], clones→1; `em_combine` with one project equals the plain VanRaden GRM. |
| `.vanraden` (`pipeline.R`) | the GRM over `need` with population frequencies **equals** the `[need, need]` submatrix of the full-population GRM; on an all-alt-major panel it equals hand-computed VanRaden while `rrBLUP::A.mat` returns all-`NaN` (the coding bug); QC on 5 accessions keeps a different marker set than QC on the 200-accession population. |
| `.group_by_panel` / `.prune_redundant` (`pipeline.R`) | 100%-marker-containment projects group together, 40% do not; identical accession sets → the richer marker build survives; ≥90%-overlapping accession sets → the project with more accessions survives; disjoint accession sets → both kept. |
| `.bridge_accessions` + em_combine (`pipeline.R`) | bridge = accessions in ≥2 panels (empty when disjoint); a panel with one needed accession is `NULL` under `need` alone but a valid GRM once bridges are kept; `build_kernel("em_combine")` is feasible *because of* the bridges and returns exactly `need × need` (bridges eliminated). |
| `build_targets` + `.blue_per_trial` (`pipeline.R`) | `raw_mean` = mean of per-trial means; `trial_center` invariant to a constant added to one trial; `blue_lm` recovers group means; outlier removal at `z_thr`; `two_stage_blup` shrinks spread; `per_trial_z` → mean 0 sd 1. |
| `predict_test` (`pipeline.R`) | conditional expectation with a clone test line → `mu + u[clone-source]`; `blend_obs_w` blends under CV0, no-op under CV00; fallback below `min_overlap` → all predictions `mean(targets)`. |
| `mask_cv` (`pipeline.R`) | CV0 keeps every row; CV00 drops exactly the focal-accession rows (CV0 set minus CV00 set = the focal rows). |

### Genome invariants *(implemented: `tests/test_config_space.R`)*

Sampling, encoding, JSON round-trips, and crossover/mutation all produce well-formed configurations. Property-based: \~5 invariants (canonical key order, valid method, param-applies⇔not-NA, JSON round-trip hash, crossover/mutation validity) asserted across \~400 randomly sampled and recombined configs — hence the `8007 passed` counter (individual assertions, not distinct cases).

### End-to-end *(implemented: `tests/test_sim_loop.R`, slow)*

Runs the whole optimizer against the synthetic world (`.sim_true`/`.sim_evaluate` in `evaluate.R`) and asserts the discovered incumbent beats the submitted seeds and that the surrogate phase improves on random search. This is the regression test for any change to the search space or optimizer; run it after such changes.

### Tiers 2–4 — planned

Specified for later; not yet implemented. Oracles already worked out:

- **Tier 2 — optimizer mechanics.** `expected_improvement` (≥0; ↑ with mean; ↑ with sd when mean\<best; →0 as sd→0); `fit_surrogate`/`predict_surrogate` (NULL below `min_obs`; monotone signal → high rank-correlation; sd floor); `aggregate_scores`/ `get_elites`/`incumbent_config` (mean excludes failures, `n` counts them; top-k; `min_reps` honored); `crossover`/`mutate_config`/`propose_candidates` (each child block matches one parent; mutation changes ≥1 block; candidates repaired); `choose_config` phases (empty store → seed; \< `n_random_init` → `random_init`; enough → `acquisition`).
- **Tier 3 — store / data helpers / reporting.** `store_eval`↔`read_evals` round-trip incl. NA score and `detail`; `.matches_trait`/`.focal_trait_parts`/ `.apply_target_domain`; `.obs_tibble`+`get_observations` on a fake response; `cached` write/read; `report.R` (`format_config`, `method_importance`, `failure_summary`, "BEATS submissions" only when incumbent \> best seed).
- **Tier 4 — network subtasks via a mock `conn`.** Inject a fake `conn` (closures for `$wizard`/`$search`/`$vcf_archived`) + a fake `trial_catalog` and test the *logic* of `select_training_trials`, `.trial_similarity`, `choose_geno_sources` (most-similar ranks first; `same_program` respects `prog_cap`; `top_k_similar` returns `k`; `best_single_project` picks the most-covering project). Needs a `tests/fixtures.R`.

------------------------------------------------------------------------

## 8. How the code touches cache / state / logs

**`cache/`** — the on-disk memoizer is `cached(settings, name, expr, max_age_days)` in `data_access.R`; every entry is `cache/<name>.rds` (or a raw `.vcf`). Touch a trial once, reuse forever.

| Cache file | Written by | Holds | Max age |
|----|----|----|----|
| `trial_catalog.rds` | `trial_catalog()` | focal-trait trials + metadata + lat/long/elev | 7 d |
| `acc_<studyid>.rds` | `get_trial_accessions()` | accession names of a trial. Also the substrate `.find_related()` and `.trial_similarity()` compute germplasm overlap from — one file per trial serves every pairwise comparison, so there is no separate neighbour cache | 30 d |
| `obs_<studyid>.rds` | `get_observations()` | all numeric observations of a study | 30 d |
| `proj_<hash>.rds` | `projects_for_accessions()` | genotyping **project** ids covering an accession set | 30 d |
| `raw_project_<id>.vcf` | `.ensure_project_vcf()` | the archived VCF (validated complete). **Transient**: deleted once its dosage is cached OR the archive is found unparseable | until dosage/unparseable cached |
| `dosage_<id>[_thin<e>]_sz<size>.rds` | `get_project_dosage()` | the **whole project's** accessions×markers dosage, parsed **once** at the densest thin `e` the `dosage_budget_bytes` budget allows (thin 1 = full markers for anything that fits). A config's requested `marker_thin` is **derived at read** by keeping every `max(1, floor(marker_thin/e))`-th marker — so a project is never re-downloaded for a different thin. One dosage file per project (older per-thin files are redundant; `prune_dosage_cache.R` removes them). **Re-densify:** on access, if *this* machine's `dosage_budget_bytes` affords a denser thin than the cached `e` (e.g. a laptop-thinned cache warmed onto a big-memory server), it re-downloads and re-parses denser once — from `stat_<id>`, cheap to detect; `dosage_redensify = FALSE` disables it. | ∞ |
| `stat_<id>.rds` | `get_project_dosage()` | `{n_samples, n_markers}` of a project, written on first extraction so a later run can compute the effective thin — and hence the right `dosage_` name — without the (deleted) VCF | ∞ |
| `unparseable_<id>.rds` | `get_project_dosage()` | **negative cache**: this archive is not a usable VCF (transposed / header-vs-data sample-count mismatch / no genotypes). Future calls skip it instead of re-downloading and re-failing every run. `{reason, when}`. **Delete to force a retry** (e.g. after the archive is fixed upstream) | ∞ |
| `calibration_lightweight.rds` | (diagnostics, when you save it) | the last calibration table | — |

**`state/`** — the single source of truth.

- `evals.sqlite` (`store.R`): one table `evals` with columns `id, config_hash, config_json, trial_id, scheme, score, n_test, status, reason, detail, seconds, ts`. Every evaluation — success **and** failure — is a row; the incumbent, surrogate training data, and "what to try next" are all recomputed from this table, so the process is killable/resumable. `status` ∈ `ok | infeasible | suspect | error`; `reason` is the machine code; `detail` is the failure funnel.
- `report.md` (`report.R::write_report`): rewritten every `checkpoint_every` iterations — incumbent, subtask-method importance, the failure log (`failure_summary`), the "⚠ Suspected bugs" section, and the running-best curve.
- `STOP`: create it (`touch state/STOP`) to halt cleanly after the current eval.

**`logs/`** — `run_optimizer.R` writes a per-iteration heartbeat to stderr; the launch command redirects it to `logs/run.out`. The startup canary result and any `CANARY ALARM` / `step error` / `FATAL` messages appear here.

------------------------------------------------------------------------

## 9. Validation & debugging tooling

### Failure classification (`R/conditions.R`)

Three deliberately distinct "didn't work" signals, so the loop knows whether to continue or stop:

- **`infeasible(code, detail, funnel, suspect)`** — this *one* (trial, config, scheme) cannot be evaluated (too little genotyped overlap, no training trials, too few markers, …). Expected. `evaluate_config_on_trial` records it (`status = "infeasible"`, `reason = code`, `detail = funnel`) and the loop samples a fresh trial + config next iteration. If the funnel shows a **cliff** (lots of data in, \~zero out), it is tagged `suspect` instead — a likely bug, broken out in the report.
- **`fatal(message, code)`** — the whole run cannot proceed (bad settings, an unimplemented method). Re-raised past `evaluate.R` so `run_optimizer()` halts instead of spinning. (Implementation note: `evaluate_config_on_trial` uses a *single* `error` handler that branches on class — re-raising a fatal from a multi-handler `tryCatch` would be swallowed by the sibling `error` handler.)
- **`sample_failed(detail)`** — no usable trial could be sampled at all. Not a config's fault; the loop tolerates a few in a row then stops (`settings$max_sample_fail`).

`failure_summary()` (`report.R`) turns the stored failures into: counts by status, dominant reasons, and the **failure rate of each subtask method** — i.e. which configurations break most often.

### The canary bug-oracle (`R/diagnostics.R`)

The "skip and continue" design means a **bug that hides data looks identical to a trial that is genuinely infeasible.** The oracle catches that, via a **calibrate-then-freeze** bootstrap against ground truth that does **not** touch our code — the five Predictathon teams' submission files.

```         
canary_anchor(settings)            -> per-trial true accession count (median of the
                                      five algorithms' CV0_Predictions.csv, read from
                                      scripts/PredictionN/submission/). No optimizer code.
calibrate_canary_trials(s, conn,    -> our pipeline's per-trial counts joined to the anchor,
                        deep=FALSE)    with a `divergent` flag. Cheap path = genotyping_projects
                                      wizard membership; deep=TRUE also extracts VCF dosage and
                                      flags a wizard-vs-dosage gap (the synonym-name signal).
print_calibration(cal)             -> the side-by-side table; investigate any `divergent` row.
canary_configs()                   -> the FROZEN coverage configs: one per focal trial, assigned
                                      so that across the canaries every method + branch-level runs.
canary_coverage()                  -> which methods/levels the configs cover (should be all).
check_canaries(s, conn)            -> run each trial under its config, every scheme. Hard
                                      CANARY ALARM on a strong trial; soft warning on a weak one
                                      (settings$canary_weak_trials). Auto-run at startup in real mode.
```

**Why this breaks the catch-22:** the oracle shares the suspect code, so it cannot *first-time-validate* it alone. The independent anchor does that: calibration runs the suspect functions **and compares each count to the anchor**; a bug rarely corrupts two independent derivations identically, so agreement is positive evidence and a divergence localizes the bug. Freeze configs only once calibration agrees (human-reviewed). After that, the frozen oracle's job is *regression* detection, which it does fine despite sharing code.

#### Reading a `check_canaries` result — bug vs data vs stale oracle

`check_canaries` is **not** a permutation sweep. It runs **one frozen config per trial**, under **both** CV schemes — 9 trials × 1 config × 2 schemes = 18 rows. The configs *differ across trials by design* (coverage): the four data-rich strong trials carry the demanding branches, the rest a light filler.

| Trial | Config it is tested under |
|---|---|
| 10673 Aurora | `em_combine` + `gblup_sommer_GE` (G×E) + `focal_plus_onehop` |
| 10675 Big6 | `top_k_similar` + `rkhs_gaussian` + `rkhs` |
| 10676 CornellMaster | `same_program` + `all_projects`(thin5) + `gblup_loo_ridge` |
| 10677 YT_Urb | `accession_overlap` + `vanRaden_single` + `gblup_rrblup` |
| 10679, 10680, + weak 10674 / 10678 / 10681 | `.canary_filler`: `best_single_project` + `vanRaden_single` + `gblup_rrblup` |

So a failing row is **config-specific**: "10673 failed" means the em_combine + sommer_GE branch didn't work *on Aurora's data*, not that the pipeline is broken. That is what makes the oracle a coverage tool — but it also means a failure can be a property of the **method × that trial's data**, not a bug.

The `status` field already encodes the bug-vs-data distinction:

| status | Meaning | Category |
|---|---|---|
| **`error`** | an *uncaught R exception* (a crash) | **code bug** |
| **`suspect`** | an infeasibility whose data funnel *cliffed* (much data in, ~0 out) | **probable data-hiding bug** |
| **`infeasible`** | the pipeline *deliberately* raised `infeasible(code)` — a data precondition (`too_few_markers`, `insufficient_geno_overlap`, `no_focal_obs`, `too_few_train_trials`, …) was not met | **data-adequacy verdict** (by design) |
| **`constant`** | predictions were produced but had *no variance* (`sd(pred)==0`) → correlation undefined → `NA` | **method × data degeneracy**, not a crash |

The `constant` case is usually a *method* limitation, not missing data: `gblup`/`direct_blup` under **CV00** masks the test accessions out of training, so `fit$u` has no entry for them and every prediction collapses to `mu` → constant. A cluster of `constant` rows that are all **CV00** is that limitation, not a bug.

Triage rule:
- `error` → **code bug** (fix it). `suspect` → **probable data-hiding bug** (localize with `diagnose_trial`).
- `infeasible` (no `suspect`) → the data *as the pipeline currently sees it* does not support that config. Could be genuine inadequacy (the T3 platform-mismatch reality), a demanding method (`em_combine` needs multi-panel data with bridge accessions a trial may lack), or an artifact of a bug that is *hiding* data (which is what the poisoning fixes addressed).
- `constant` (esp. CV00) → **method degeneracy**, expected for some method × scheme pairs.

**Two caveats that make a stale result the third category, beyond bug vs data:**
1. **The oracle is frozen.** Its promise ("these configs work on these trials") was set at calibration time. If the T3 data changed — or was *hidden* by a since-fixed bug (poisoned caches, 401s) — the frozen configs can report `infeasible`/`constant` with **no code bug and no true data change**: the freeze has simply drifted. Re-calibrate (`calibrate_canary_trials` → review → re-freeze `canary_configs`) after any change that alters what data the pipeline sees.
2. **A result taken before a data-layer fix cannot be trusted as a data verdict.** E.g. a `check_canaries` run predating the cache-poisoning / auth-retry fixes may show `infeasible`/`constant` for trials whose data was being hidden — re-run on the fixed code before concluding "the data is inadequate."

### `sweep_rich_trials(settings, conn)` — the CODE-correctness oracle (`R/diagnostics.R`)

The complement to `check_canaries`, built to escape its central ambiguity. `check_canaries` runs **one config per trial across nine trials** — a *coverage* test that conflates a code bug with data-adequacy (a demanding config can legitimately fail on a data-poor trial). `sweep_rich_trials` instead **varies one method/branch at a time from a robust baseline** and runs each variant on a couple of **data-rich** trials (default Big6 `10675` + YT_Urb `10677`, via `settings$oracle_trials`). The point: on trials rich enough to satisfy every method's data prerequisites, **feasibility is a property of the code, not the data** — so a branch that can't produce a prediction on *either* rich trial is a genuine bug signal, not a data verdict.

- **22 variants** (`.oracle_variants()`) — one per subtask method plus each behaviour-changing branch level; a test asserts they cover every method in `SUBTASKS`.
- **Verdict per variant:** an `error` (crash) anywhere is always `BUG`; otherwise the variant is `ok` if it produced `ok` on ≥1 **non-degenerate** (trial, scheme), and `SUSPECT` if it never did.
- **Accuracy is reported.** The per-cell table now includes the `score` (the Pearson *r* from `score_predictions`), and the returned `$verdict` carries `best_score` (best *r* over a variant's ok cells) — so you see how *well* each branch predicts, not just that it ran. The full per-cell `score`/`n_test` are in the returned `$rows`.
- **Known-degenerate cells are excluded** (`.oracle_degenerate`): `direct_blup` under **CV00** masks the test lines out of training, so predictions collapse to the mean → `constant` by construction, not a bug. (cond_expectation predicts them from the kernel, so it is *not* degenerate.)
- **Run it:** `sweep_rich_trials(s, conn)`. `SWEEP CLEAN` = every branch predicted on rich data; a `SWEEP ALARM` lists the suspect variants — cross-check one with `diagnose_trial()` on the same trial+config. Costs ~one download of the two trials' genotypes (cached thereafter) plus a fit per variant × 2 trials × 2 schemes.
- **Re-run one branch:** `sweep_rich_trials(s, conn, only = "em_combine")` runs just the variants whose label contains the string (e.g. `only = c("kernel=", "model=")`), so after a fix you re-check one branch without the whole sweep.

This is the trustworthy, achievable half of validation (the full all-configs × all-trials matrix is neither): a green sweep says the code paths work; it does **not** claim every config suits every trial (nor should it).

### `diagnose_trial(id, settings, conn)` (`R/diagnostics.R`)

Replays one trial and prints the data funnel stage by stage **next to an independent re-derivation of the raw counts from T3** — including, per genotyping project, whether the dosage matrix rownames actually intersect the trial's accession names. A line like `overlap with accessions = 0` on a non-empty VCF is the decisive synonym/name-mismatch signal. Reach for it whenever `subtaskC` (§L8) or the `flow` funnel (§L11) looks wrong.

### Lessons baked in (don't re-learn the hard way)

- **Cache poisoning + big VCFs.** A partial/raced VCF download once cached a tiny dosage matrix *forever* (`get_project_dosage` used `max_age_days = Inf`). Now `.ensure_project_vcf()` validates completeness (`.vcf_complete`) and re-downloads truncated files. `get_project_dosage()` downloads and parses each project **exactly once**, at the densest thin the budget allows, and reuses that one cache for every trial and every requested `marker_thin` by subsetting rows (samples) and columns (markers) at read; it **deletes the raw VCF once the dosage is cached**. If genotype counts ever look low, clear `cache/dosage/` (regenerable).
- **Interactive prompt.** `conn$vcf_archived()` can prompt to pick a file and hang a non-interactive run — always launch real-mode scripts with `</dev/null`.
- **The VCF parser is streaming, and some archives are broken.** `.vcf_to_dosage()` streams the VCF in chunks (base R, no `vcfR`) so memory is bounded by one chunk plus the thinned result — a 7.5M-marker panel no longer OOMs. It (a) decompresses gzip/BGZF transparently, (b) **skips a malformed variant line** rather than failing the whole file (a real archive had one 3-field record), and (c) **rejects** a transposed / non-VCF layout. Real T3 archives seen that are genuinely unusable — a transposed export, and one whose `#CHROM` header declares more samples than the data rows carry — are recorded in `cache/unparseable/unparseable_<id>.rds` and skipped forever after (delete that file to retry). A project too large to hold densely is thinned to fit `settings$dosage_budget_bytes` at parse time, the thin `e` recorded in the `dosage_` filename; that budget is the sole control on cached marker density. If a project you expect silently yields no genotypes, look for its `unparseable_` marker and read the reason.
- **A timing-out download is not re-stormed.** A transient VCF-download failure (T3 choking on a big archive, e.g. project 8217) is tracked **per session** in RAM by project id: each subsequent attempt tries less hard (fewer in-call retries, shorter backoff) and after `settings$vcf_max_download_attempts` (default 3) the project is **skipped for the rest of the run** instead of re-downloading on every covering trial × scheme. This is *not* an `unparseable_` (structural) verdict — it is session-only and resets on a new run, so a project that was merely down is retried fresh next time. A successful download clears the counter. If a project you expect is being skipped, it hit its download-attempt cap this session — start a new session (or raise `vcf_max_download_attempts`) once T3 is healthy.

------------------------------------------------------------------------

## 10. Extending the search space

To add a literature method:

1.  Add the method name (and any new parameters, tagged with the methods they apply to) to the relevant subtask in `R/config_space.R::SUBTASKS`.
2.  Add a matching dispatcher branch in `R/pipeline.R` (the `select_training_trials` / `build_targets` / `choose_geno_sources` / `build_kernel` / `train_model` / `predict_test` switch for that subtask). An unhandled method must `fatal()`.
3.  Re-run `Rscript tests/run_all.R` (genome invariants + subtasks) and, after a search-space change, `tests/test_sim_loop.R`.
4.  Re-run `canary_coverage()` so the new method is exercised by a canary (extend a `canary_configs()` entry if needed).
5.  If you added a function worth stepping through, add its name to the appropriate group in `R/evaluation.R::EVAL_GROUPS` so `arm_evaluation()` reaches it.

Nothing else in the optimizer needs to change — the surrogate encoding, crossover, mutation, store, and report all derive from `SUBTASKS`.
