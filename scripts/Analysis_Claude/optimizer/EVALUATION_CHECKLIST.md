# Optimizer evaluation checklist

Tick these off over time. Each line names the `arm_evaluation()` group and the full step-by-step walkthrough lives in **`EVALUATION.md`** at the matching level (L1…L12). Order is fast/offline → slow/online on purpose: do it top to bottom so a cheap bug surfaces before an expensive one. `disarm_evaluation()` at the end of each level, and **always before launching a background run**.

Bootstrap once (see EVALUATION.md §3), then:

## Offline — `simulate = TRUE`, no network, do these first

- [x] **L1 `genome`** — `sample_config()` / `crossover()` / `mutate_config()` produce well-formed configs; `config_hash()` stable; inapplicable params are `NA`, applicable ones are not; a mutation changes ≥1 block.
- [x] **L2 `engine`** — `choose_config()` moves through phases seed → `random_init` → `acquisition` at the right counts; `expected_improvement()` ≥ 0; on the sim world the surrogate phase beats random search; incumbent honours `incumbent_min_reps`.
- [x] **L3 `store`** — `store_eval` ↔ `read_evals` round-trip (incl. `NA` score & `detail`); `config_from_json(config_to_json(cfg))` == `cfg` including `NA` params; `write_report()` renders `state/report.md`.
- [x] **L4 `scoring`** — `score_predictions(x, a·x+b)` = +1 (a\>0); \<5 overlap or constant → `NA` with a reason; `mask_cv` under CV00 drops **exactly** the focal accessions CV0 keeps; no leakage (score not ≈1).
- [x] **Test suite** — `Rscript tests/run_all.R` → `2/2 test files passed`; counts match (`config_space` 8007, `subtasks` 196). `run_all.R --all` → sim-loop `PASS`.

## Live BrAPI — `simulate = FALSE` — smallest / cheapest first

- [x] **L5 `data`** — on the **cached** trial `10676` (no download): descriptor fills; `peek(obs)` shows finite `value` and **`rep`/`block` not all `NA`**; `n_acc` matches the T3 web UI.
- [x] **L6a `subtaskA` / `accession_overlap`** — a **HIGH** `primary_min` (few trials) before a low one; result respects the threshold; not unexpectedly empty; the secondary tier is additive (primary ⊆ primary+secondary) and actually adds something.
- [x] **L6b `subtaskA` / `top_k_similar`** — returns ≤ `k`; raising `k` is nested; the `genomic` ranking matches an independent shared-accession count. **Check `is.finite(trial$lat)` before trusting `environmental`/`both`** — without coordinates every score is 0 and the `k` returned are arbitrary.
- [x] **L6c `subtaskA` / `same_program`** — every id is in the focal program, focal excluded, length ≤ `prog_cap`. If it comes back **empty**, check the program's size first: the `> prog_cap` same-location-or-year fallback can keep nothing.
- [x] **L7 `subtaskB`** — `peek(train_obs)` no column `ALL NA`; `peek(targets)` distinct ≫ 1 and a plausible yield range (not degenerate).
- [x] **L8 `subtaskC`** — `peek(dl, accessions=…)` every dosage matrix has **rowname overlap \> 0** (synonym check); `n_projects > 0 ⇒ length(dl) > 0`. Cross-check `diagnose_trial("10676", s, conn)`.
- [x] **L9 `subtaskD`** — `peek(K)` square, **symmetric**, sane diagonal; clones give identical rows; ≥50 markers survived QC.
- [x] **L10 `subtaskE` + `subtaskF`** — `peek(pred)` not constant (when overlap adequate); CV0 blend on, CV00 blend is a no-op.
- [x] **L11 `flow`** — `run_pipeline()` funnel has **no cliff** (many in, \~0 out); `score_predictions(po$pred, po$obs)` plausible (not ≈1, not `NA`).
- [ ] **L12 `diagnostics`** — `calibrate_canary_trials()` agrees with `canary_anchor()` (no `divergent` rows — the light, trustworthy plumbing check); `sweep_rich_trials()` is `SWEEP CLEAN` (the code oracle — every branch predicts on rich trials); `check_canaries()` green on strong trials (coverage oracle; weigh `infeasible`/`constant` per §9's bug-vs-data table, soft-warn OK on weak trials `10674/10678/10681`).

## Whole system

- [ ] **Short real run** — `run_optimizer(modifyList(optimizer_settings(), list(simulate = FALSE, max_iters = 10)), conn)` completes, checkpoints `state/report.md`, and the incumbent ≥ the best seed.
- [ ] **Clean shutdown** — after evaluating, `disarm_evaluation()` restores a debug-free state so a background run will not block on the debugger.

------------------------------------------------------------------------

*Notes / anomalies noticed (the "a little funny" observations — the point of the exercise):*

- 

- 

- 
