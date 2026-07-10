# Optimizer evaluation checklist

Tick these off over time. Each line names the `arm_evaluation()` group and the
full step-by-step walkthrough lives in **`EVALUATION.md`** at the matching level
(L1…L12). Order is fast/offline → slow/online on purpose: do it top to bottom so a
cheap bug surfaces before an expensive one. `disarm_evaluation()` at the end of
each level, and **always before launching a background run**.

Bootstrap once (see EVALUATION.md §3), then:

## Offline — `simulate = TRUE`, no network, do these first

- [ ] **L1 `genome`** — `sample_config()` / `crossover()` / `mutate_config()` produce well-formed configs; `config_hash()` stable; inapplicable params are `NA`, applicable ones are not; a mutation changes ≥1 block.
- [ ] **L2 `engine`** — `choose_config()` moves through phases seed → `random_init` → `acquisition` at the right counts; `expected_improvement()` ≥ 0; on the sim world the surrogate phase beats random search; incumbent honours `incumbent_min_reps`.
- [ ] **L3 `store`** — `store_eval` ↔ `read_evals` round-trip (incl. `NA` score & `detail`); `config_from_json(config_to_json(cfg))` == `cfg` including `NA` params; `write_report()` renders `state/report.md`.
- [ ] **L4 `scoring`** — `score_predictions(x, a·x+b)` = +1 (a>0); <5 overlap or constant → `NA` with a reason; `mask_cv` under CV00 drops **exactly** the focal accessions CV0 keeps; no leakage (score not ≈1).
- [ ] **Test suite** — `Rscript tests/run_all.R` → `2/2 test files passed`; counts match (`config_space` 8007, `subtasks` 45). `run_all.R --all` → sim-loop `PASS`.

## Live BrAPI — `simulate = FALSE` — smallest / cheapest first

- [ ] **L5 `data`** — on the **cached** trial `10676` (no download): descriptor fills; `peek(obs)` shows finite `value` and **`rep`/`block` not all `NA`**; `n_acc` matches the T3 web UI.
- [ ] **L6 `subtaskA`** — `select_training_trials` with a **HIGH** `accession_overlap` `primary_min` (few trials) before a low one; result respects the threshold; not unexpectedly empty.
- [ ] **L7 `subtaskB`** — `peek(train_obs)` no column `ALL NA`; `peek(targets)` distinct ≫ 1 and a plausible yield range (not degenerate).
- [ ] **L8 `subtaskC`** — `peek(dl, accessions=…)` every dosage matrix has **rowname overlap > 0** (synonym check); `n_projects > 0 ⇒ length(dl) > 0`. Cross-check `diagnose_trial("10676", s, conn)`.
- [ ] **L9 `subtaskD`** — `peek(K)` square, **symmetric**, sane diagonal; clones give identical rows; ≥50 markers survived QC.
- [ ] **L10 `subtaskE` + `subtaskF`** — `peek(pred)` not constant (when overlap adequate); CV0 blend on, CV00 blend is a no-op.
- [ ] **L11 `flow`** — `run_pipeline()` funnel has **no cliff** (many in, ~0 out); `score_predictions(po$pred, po$obs)` plausible (not ≈1, not `NA`).
- [ ] **L12 `diagnostics`** — `calibrate_canary_trials()` agrees with `canary_anchor()` (no `divergent` rows); `check_canaries()` green on strong trials (soft-warn OK on weak trials `10674/10678/10681`).

## Whole system

- [ ] **Short real run** — `run_optimizer(modifyList(optimizer_settings(), list(simulate = FALSE, max_iters = 10)), conn)` completes, checkpoints `state/report.md`, and the incumbent ≥ the best seed.
- [ ] **Clean shutdown** — after evaluating, `disarm_evaluation()` restores a debug-free state so a background run will not block on the debugger.

---

_Notes / anomalies noticed (the "a little funny" observations — the point of the exercise):_

-
-
-
