# Lessons — traps we fell into, and what stops us falling in again

Every entry here is something that **looked correct, ran without error, and produced a plausible
wrong answer**. They are recorded away from the code because a fresh reader does not need our
history to understand what the code does — but *we* need it, because most of these are re-enterable
by a well-meaning simplification.

The code cites these by number (`see docs/LESSONS.md #11`). Numbers are stable: append, never renumber.

**This is not the debugging guide.** When something looks wrong *right now*, go to `docs/EVALUATION.md`
§6 (subtle-bug catalogue) and §9 (tooling). This file is why those checks exist.

| | |
|---|---|
| §1 The T3 data layer | #1–#10, #27 |
| §2 Genomic methods | #11–#14 |
| §3 Search space and design | #15–#17, #33 |
| §4 Testing and validation | #18–#20 |
| §5 Numerical robustness | #21–#22 |
| §6 Running it: resources, concurrency, configuration | #23–#26, #28–#32, #34–#45 |

---

## §1. The T3 data layer

Getting clean inputs out of T3/Wheat over BrAPI is a substantial part of the work: responses are
inconsistent, and several plausible-looking calls return wrong-but-non-erroring results. Fixes are
in `R/data_access.R` unless noted.

### 1. Trait → trial filtering needs the numeric variable id, not the name

**Trap.** Passing the human trait string (`"Grain yield - kg/ha|CO_321:0001218"`) to the breeder
wizard's `traits` filter.
**Symptom.** *Every* trial comes back, with only a "no matching trait" warning — the wizard's trait
dimension lists only derived `COMP:` traits. Silently, the optimizer would be sampling focal trials
that never measured the trait.
**Now.** The studies search is filtered by `observationVariableDbIds` (`settings$focal_trait_db_id`,
grain yield = `84527`), narrowing ~9030 wheat trials to ~7561. `focal_trait` (the name) is still
used for matching *within* an observation set.
**Lives in.** `trial_catalog()`, `settings$focal_trait_db_id`.

### 2. Trial metadata carries no coordinates and no year

**Trap.** Assuming `get_all_trial_meta_data` returns what the environmental paths need.
**Symptom.** `trial$lat` not finite → the whole environmental-similarity branch is skipped, every
score stays 0, and `top_k_similar` returns an *arbitrary* k trials that look like a real answer.
**Now.** `year` is derived from `start_date`; `latitude`/`longitude`/`elevation` are joined per
location via `T3BrapiHelpers::get_lat_long_elev_from_location_vec` and cached with the catalogue.
This is what makes the environmental-similarity and G×E paths possible at all.
**Lives in.** `trial_catalog()`, `.trial_similarity()`.

### 3. snake_case and camelCase are not interchangeable

**Trap.** Mixing the two column conventions in one pipeline.
**Symptom.** "object not found" inside `dplyr::filter` — or worse, a silent empty filter.
**Now.** All `T3BrapiHelpers` outputs are `janitor::clean_names()`-ed (`study_db_id`,
`program_name`, …); raw `conn$search` / `conn$wizard` responses keep camelCase. Check which side of
that line a value came from before joining.

### 4. Genotyping protocol id ≠ project id

**Trap.** Feeding `get_geno_protocol_from_germ_vec`'s output to `conn$vcf_archived()`.
**Symptom.** No downloadable archive; a trial looks ungenotyped when it is not.
**Now.** Project ids come from the `genotyping_projects` wizard (`projects_for_accessions()`), with
accessions batched.

### 5. `/observations` hides its payload

**Trap.** Reading `resp$data`.
**Symptom.** `resp$data` is search metadata, not records — an empty or nonsensical phenotype set.
**Now.** Records are read from `resp$combined_data` (the auto-paginated list), each a nested list.
Note that `rep`/`block` are *not* returned (they would need an `observationUnits` query), so BLUE
estimation degrades gracefully to per-germplasm means — which is legitimate, but means an all-`NA`
`rep` column is expected rather than alarming.

### 6. VCF sample names may be synonyms

**Trap.** Matching focal accessions to VCF sample names directly.
**Symptom.** A non-empty dosage matrix with *zero* rowname overlap — an accession is genotyped, yet
invisible. Reads downstream as ordinary infeasibility.
**Now.** Sample names are canonicalized to primary germplasm names before matching. **Only partially
handled**: a residual tail remains (~22 of 139 on one trial), and `calibrate_canary_trials()`
measures its size per trial. This is the one open item in this section.

### 7. Cache poisoning from a partial download

**Trap.** Caching a parsed artefact with `max_age_days = Inf` without validating its source.
**Symptom.** A truncated/raced VCF download produced a tiny dosage matrix that was then cached
*forever*, silently shrinking every GRM built from it. Observed: one project's focal overlap
collapsed from 115 accessions to 6.
**Now.** `.ensure_project_vcf()` validates completeness (`.vcf_complete`) and re-downloads a
truncated file; the dosage cache key includes the VCF byte size **and** a hash of the requested
samples, so a re-downloaded VCF or a different sample set can never return a stale subset. If
genotype counts ever look low, clear `cache/dosage/` — it is regenerable.

### 8. `vcf_archived()` can prompt, and a prompt hangs a batch run

**Trap.** Running real-mode code non-interactively without redirecting stdin.
**Symptom.** The job hangs indefinitely with no error, waiting for a file choice.
**Now.** Always launch real-mode scripts with `</dev/null`.

### 9. Some archives are not usable VCFs, and big ones must stream

**Trap.** Loading a whole VCF, and assuming every archive is well-formed.
**Symptom.** OOM on a 7.5M-marker panel; or one malformed record killing an entire file.
**Now.** `.vcf_to_dosage()` streams in chunks (base R, no `vcfR`), so memory is bounded by one chunk
plus the result. It decompresses gzip/BGZF transparently, **skips** a malformed variant line rather
than failing the file (a real archive had a 3-field record), and **rejects** a transposed/non-VCF
layout. Genuinely unusable archives — a transposed export, and one whose `#CHROM` header declares
more samples than the data rows carry — are recorded in `cache/unparseable/unparseable_<id>.rds` and
skipped thereafter. Delete that file to retry. If a project you expect yields no genotypes, look for
its `unparseable_` marker and read the reason.
**Lives in.** `R/genotypes.R::.vcf_to_dosage` / `.vcf_header` / `.vcf_stat`, `.neg_cache_mark()`.

### 10. A timing-out download must not be re-stormed

**Trap.** Retrying a failed download with full effort on every trial that needs it.
**Symptom.** One sick archive (e.g. project 8217) consumes an entire run, re-downloading on every
covering trial × scheme.
**Now.** Failures are tracked **per session** in RAM by project id: each attempt tries less hard
(fewer retries, shorter backoff), and after `settings$vcf_max_download_attempts` (default 3) the
project is skipped for the rest of the run. This is deliberately *not* an `unparseable_` verdict —
it is session-only and resets next run, so a project that was merely down is retried fresh. A
success clears the counter.
**Lives in.** `R/genotypes.R::.vcf_download_plan` / `.ensure_project_vcf`.

### 27. A rejected password logs in "successfully"

**Trap.** Calling `conn$login()` and assuming a wrong password produces an error.
**Symptom.** BrAPI's `login()` assigns `resp$content$access_token` unconditionally, so a
REJECTED password leaves `auth_token` NULL and **returns normally** — with only an "Incorrect
Password" *warning* in the output. Every later call then goes out anonymous and 401s, and the
run reports data-shaped failures: no descriptor, empty searches, trials that look ungenotyped.
Cached counts stay healthy, so the store looks fine while nothing new arrives. A one-line
credentials problem presents as a data problem, for as long as you let it run.
**Now.** `t3_login()` verifies a token was actually issued and raises `t3_bad_credentials`,
which `.brapi_try()` treats as fail-fast rather than transient — retrying would only re-read
the same environment. Missing credentials raise the sibling `t3_missing_credentials`.
**Generalise.** For any auth call, assert on the *artefact* (a token, a session id), never on
the call returning without error.
**Lives in.** `R/data_access.R::t3_login`, `.brapi_try`.

---

## §2. Genomic methods

### 11. `rrBLUP::A.mat()` silently destroys a {0,1,2} dosage matrix

**Trap.** Reaching for the library GRM function. It is the obvious, idiomatic choice.
**Symptom.** No error. `A.mat` documents `{-1,0,1}` coding; fed `{0,1,2}`, its internal
`freq = (colMeans(X)+1)/2` equals *p* + 0.5, so `MAF = pmin(freq, 1-freq)` goes **negative** for
every marker whose alternate allele is the major one, and `min.MAF` drops them — roughly half the
panel. The resulting GRM's off-diagonals correlate only ~0.7 with the correct one. Centering happens
to survive; the scaling does not.
**Now.** `.vanraden()` computes the GRM directly. The `{0,1,2}` assumption did not disappear — it
*moved* into `p`, so `.vanraden()` guards it: a negative entry means someone changed the dosage
encoding, and it fails loudly rather than reproducing the same class of silent error. A companion
oracle test asserts `A.mat` returns all-`NaN` on an all-alt-major panel where `.vanraden` is correct.
**Lives in.** `R/pipeline.R::.vanraden`, `tests/test_subtasks.R`.

### 12. Estimate marker QC and allele frequencies on the population, not the request

**Trap.** Running QC over just the accessions a trial needs — it is the smaller matrix, so it looks
like an optimization.
**Symptom.** The surviving marker set changes with *who asked*: 10 needed accessions kept ~6.9k
markers where the panel population kept ~9.0k, and allele frequencies estimated from a handful of
lines are unstable. Two trials then get incomparable GRMs from the same panel.
**Now.** Subtask C returns each panel's **whole population**; subtask D runs QC and estimates
frequencies on that, then forms the GRM over the needed accessions. Since `K_ij` depends on other
accessions only through *p*, this is exactly the `[need, need]` submatrix of the population GRM, at
`O(n_need²·m)` rather than `O(n_pop²·m)` — correct *and* cheaper.
**Lives in.** `R/pipeline.R::choose_geno_sources`, `build_kernel`, `.qc_markers`.

### 13. Bridge accessions are scaffolding — drop them too early and the combine collapses

**Trap.** Building each per-panel GRM over only the accessions the trial needs.
**Symptom.** `em_combine`'s cross-panel blocks come back ~0 even when the panels *do* share
accessions, so it silently degenerates to a block-diagonal matrix — the entire value of the method,
gone, with no error.
**Now.** Per-panel GRMs retain the **bridge** accessions (genotyped in >1 protocol group) so
`covariance_combiner()` has shared rows to stitch through; the bridges are dropped from the combined
result, which is subset back to `need`.
**Lives in.** `R/pipeline.R::build_kernel`, `.bridge_accessions`.

### 14. Protocol ids do not identify a marker panel

**Trap.** Grouping genotyping projects by protocol id.
**Symptom.** A `V2` and a `v2.1` protocol are the *same* assay scored against a different reference
genome: different ids, near-identical markers. Grouping by id splits a panel in two and merges
nothing.
**Now.** Projects are grouped by **marker overlap** (`settings$merge_containment`, 0.95). On the T3
projects seen so far the classes separate cleanly: same-panel pairs share 100% of the smaller panel,
everything else ≤75%.
**Lives in.** `R/pipeline.R::.group_by_panel`.

---

## §3. Search space and design

### 15. A declared knob is not a wired knob

**Trap.** `SUBTASKS` (`R/config_space.R`) and the dispatchers (`R/pipeline.R`) are two files, and
nothing enforced that an entry in the first had a branch in the second.
**Symptom.** For a long time three did not. `geno_select.min_bridge` was sampled and never read, so
`focal_plus_onehop` was byte-for-byte `all_projects`. `model.lambda_select` was sampled and never
read. `model.method = "gblup_loo_ridge"` fell through to the plain rrBLUP backbone. Nothing failed:
the optimizer searched over distinctions that did not exist, wasting evaluations on duplicates in
disguise — and the SIMULATE world *scored them as if they were real*, so the offline tests stayed
green throughout.
**Now.** All three are implemented. The sweep-coverage test asserts every method in `SUBTASKS` is
exercised by an oracle variant, which would have caught the dead *method* (though not the dead
*parameters*). When adding to the search space, follow `docs/EVALUATION.md` §10 in full. And when a
method looks strangely indistinguishable from another in `method_importance()`, check that it *is*
distinguishable in the code before believing the result.

### 16. A knob can be *read* and still be inert

**Trap.** `geno_select.marker_thin` was genuinely wired — read by `choose_geno_sources`, applied by
`get_project_dosage` — and was still removed.
**Symptom.** A lower layer had already made the choice. Each project is parsed once at the densest
thin the memory budget allows, and a request could only *subset* that cache, so every level returned
an identical matrix whenever `dosage_budget_bytes` had thinned harder than the request — true for
exactly the huge panels thinning was meant to help. It also made a configuration's meaning depend on
which machine ran it.
**Now.** Marker density is a settings-level property (`dosage_budget_bytes`), not a search
parameter. The tell to watch for: a parameter whose effect depends on a **settings** value rather
than on the configuration. Search only over what the pipeline can vary on its own.

### 17. A convenient helper can be quadratically expensive

**Trap.** `T3BrapiHelpers::find_other_studies_evaluating_same_germplasm()` is the obvious way to find
neighbouring trials.
**Symptom.** It spends one wizard query **per germplasm** — 223 queries to answer for a single
trial — and answers for one trial only.
**Now.** `.find_related()` issues one wizard query for candidate trials and computes overlaps as set
intersections of the `acc_<sid>` cache, which is shared across every trial and with
`.trial_similarity()`, so the cost decays to nothing as the cache warms. Two deliberate differences:
overlap is counted in germplasm **name** space (what every downstream join uses), and candidates are
restricted to the focal-trait catalogue — a trial that never measured the trait is not a training
trial.
**Lives in.** `R/pipeline.R::.find_related`.

### 33. The search can churn: turnover among candidates with no improvement and no resolution

**Trap.** Reading a plateau in the leader's score plus constant turnover in the contender list as
"the search is still exploring", and waiting for replication to settle the field.
**Symptom.** On run `9a2664962c01` the leader's estimate rose from 0.324 to 0.342 by evaluation
703 and then sat between 0.336 and 0.339 for the remaining 5,533 evaluations, while the leader's
identity changed at 19 of 25 checkpoints (19 distinct leaders, median top-8 Jaccard 0.33) and the
leader's usable replication was often 1–3. Nothing was converging and nothing was being decided.
Run `80cfbd7c27e9`, same code and build, did not churn (5 leaders, Jaccard 0.89) — so this is a
regime, not a constant, and a run cannot be judged by another run's behaviour.
**The five mechanisms, each measured.** (1) `sd_trial` exceeds `sd_config` in six of the seven
stored runs. (2) `se` can never exceed `sd_config`, and replication shrinks it only in proportion
to `n/n*` where `n* = (sd_resid/sqrt(w_bar))^2 / sd_config^2`; on the two 24-trial targets `n*` is
24.9 and 26.6, i.e. larger than the universe, so the best precision a full budget can reach there
is `se = 0.7 sd_config` — the posteriors stay overlapped whatever the budget.
(3) The contender cap is the top 8 by UCB, with no memory: a configuration pushed out is in no
backlog tier, is not re-proposed by acquisition, and is frozen — and because the extra evaluation
pulls a top-ranked BLUP back toward the grand mean, replication *removes* configurations from the
set rather than accumulating evidence on them (the previous checkpoint's leader sat at median rank 19 at the
next). (4) The EI uncertainty floor bound for 121 of 121 candidates, which equalises their
uncertainty and makes EI a strictly increasing function of the predicted mean — plain greedy
exploitation — and since the forest predicts lower for fresh random draws than for elite
recombinations, exploration stops; 97% of that run's 1,043 configurations then shared five
identical continuous parameters and 60 crossover draws de-duplicated to 2. (5) That homogeneous
population drives `sd_config` down (0.0263 → 0.0078 over the run; exactly 0.0000 when the second
half is fitted alone), which strengthens the shrinkage in (3) and the greediness in (4).
**Now.** Nothing is fixed. What exists is the diagnosis: `docs/SEARCH_ALGORITHM.md` defines churn
as three required series (turnover, direction, resolution), gives `n*` against the universe size
as the cheap test of how much precision a target allows at all, and documents each mechanism against
its code. `BACKGROUND.md` holds the candidate remedies (shared trial panel, empirical-Bayes
blending, racing in place of the fixed cap), none built.
**Generalise.** A quantity that is both the prior and the yardstick will dissolve the decision
rule it is supposed to enforce. `sd_config` is the shrinkage prior *and* the unit the contender
threshold is measured in, so a search that concentrates its population automatically shrinks the
scale on which "still in play" is judged. A threshold worth trusting is denominated in a fixed,
externally chosen difference — the `Delta` in predictive ability actually worth resolving — not in
a variance the search itself moves.
**Lives in.** `R/optimizer.R::.contenders` / `.replication_backlog` / `.blup_scores`,
`R/surrogate.R::predict_surrogate`.

---

## §4. Testing and validation

### 18. A pass/fail assertion built on a noisy statistic tests the RNG

**Trap.** `test_sim_loop.R` asserted that the mean score of the last 50 evaluations exceeded the mean
of the first 25 — an intuitive reading of "the search improves."
**Symptom.** Every evaluation lands on a *different* synthetic trial, and trial difficulty dominates
the config's contribution, so the comparison held in only ~7 seeds out of 8. Removing an unrelated
search-space parameter reshuffled the random stream and flipped it — a green suite turning red with
no regression, which is the same failure as a red suite staying green.
**Now.** Two questions are kept apart. *Does the algorithm work?* is a correctness property, asserted
in a **deterministic regime** (`sim_noise_sd = 0`, `sim_fixed_trial = TRUE`) where a correct
implementation must beat random search — it does, 0.606 vs 0.553, in every seed. *How many
evaluations does a given noise level need?* is an empirical property of the objective, measured and
recorded, asserted nowhere. The general rule: assert correctness where correctness is decidable, and
measure performance where it is merely quantifiable.

### 19. Trial heterogeneity, not measurement noise, is what defeats the search

**Trap.** Diagnosing "the optimizer is only nominally better than random" as observation noise, and
reaching for more precise evaluations.
**Symptom.** With observation noise switched **off** and trials still varying, the search is *still*
no better than random (0.513 vs 0.525). The variance that matters comes from the objective's own
definition — a mean over a heterogeneous trial population, estimated one trial at a time — not from
measurement.
**The measurement, on the real store.** `sd_trial` **0.078** against `sd_config` **0.036** — trials
differ more than twice as much as the thing being optimized, and idiosyncratically: program, year,
location and test-set size explain essentially none of it. That is why the surrogate blocks on
`trial_id` itself rather than on trial descriptors, and why blocking needs replication to be safe —
at one observation per trial an `rpart` split on `trial_id` fits pure noise to an in-sample R² of
0.98.
**Now.** Fixed, in two halves. `trial_replication` makes a trial carry several configurations, so
the trial effect is estimable and `aggregate_scores()` removes it with random-effects BLUPs.
`config_replication` (0.8.5) makes a configuration carry several trials, so its mean is an average
over the population it is supposed to generalize to. The two cooperate: the trial backlog supplies
exactly the under-used trials a replication step needs. Both levels **ramp** (0.8.6) — trials on a
`sqrt(n_scored)` schedule, configurations by contention, so the budget concentrates on whichever
configs are still in question. See `docs/DESIGN.md` §3 and `dev/README.md` (Open threads).

**Read the variance components carefully.** `.blup_scores()` weights by `n_test - 3`, so lmer's
`sd_resid` is the residual sd at *unit weight* while `sd_config` and `sd_trial` are on the score
scale. Comparing the three as printed makes the residual look catastrophic; divide it by
`sqrt(median(n_test) - 3)` first. `report.md` now prints both.

### 20. Your analysis scripts deserve the same suspicion as the pipeline

**Trap.** Trusting a throwaway measurement script because the code under test is the "real" code.
**Symptom.** A comparison of two optimizer arms reported differences of *exactly* zero across nine
cells. The cause was a duplicated `n_random_init` key in a `list()` passed to `modifyList`, which
takes the first match — so both arms ran the same configuration. The tell was not an error but an
implausibly clean result.
**Now.** The comparison script asserts its own arms differ before running. The general rule is the
one behind every entry here: **trust agreement between two independent derivations; distrust a
single number** — including one your own analysis produced.

---

## §5. Numerical robustness

### 21. An EM combiner needs positive-definite inputs, and the ridge was added too late

**Trap.** `build_kernel()` built a GRM per panel, scaled each to a unit mean diagonal, handed them
to `covariance_combiner()`, and added `ridge` to the diagonal of the **combined** result. Reading
top to bottom it looks regularized.
**Symptom.** The EM inverts those partials, and a raw GRM is easily singular — identical rows from
clones or from two names for one line, or fewer surviving markers than kept accessions. Three
evaluations died with `Lapack routine dgesv: system is exactly singular` and reciprocal condition
numbers around 1e-16, one of them after 13.4 hours of work. `status = "error"`, so it read as a
code bug, which it was.
**Now.** Each partial is ridged (with a floor, so `ridge = 0` is still safe) *before* the combine;
a combine that fails anyway raises `infeasible` carrying every partial's shape and rank, rather
than an uncaught exception. Rank-deficient partials are reported once per session with their
duplicate-row count — a non-zero count is evidence about the synonym tail (#6), so the
regularization does not hide the data problem it papers over.
**See also.** `dev/EM_COMBINE_COMPARISON.md` — the sibling `Brapi_pipeline_for_selection` also
regularizes *inside* the E-step (`psi_aa + diag(1e-5)`), which `T3BrapiHelpers` does not, so the
one place we cannot ridge from outside is exactly where it does. That document also covers the
degrees-of-freedom difference (this pipeline weights panels by accession count; the sibling by an
LD-aware effective n) and pedigree bridging, which bears directly on #13.
**Lives in.** `R/pipeline.R::build_kernel` (the `em_combine` branch).

### 22. A method that silently changes meaning is worse than one that fails

**Trap.** `predict_test()`'s `direct_blup` took `fit$u[test_in]` and replaced `NA` with 0. Sensible
defensive coding.
**Symptom.** The rrBLUP backbone is fitted on the training lines alone, so under CV00 — where the
focal lines are masked out — it has no effect for them at all. Every prediction became `mu`: a
constant vector, an `NA` score, and no error anywhere. Nine of sixteen evaluations in one window
ended this way, some after 100+ minutes. The optimizer was buying guaranteed-unscoreable
evaluations at an hour each, and four of the five seed configurations were affected, so the CV00
baseline barely existed.
**Now.** When the fit carries no effect for the focal lines, prediction goes through the kernel —
which is what all five submitted algorithms do for masked lines (three by explicit conditional
expectation, two by carrying the test lines as levels in the fitted model) — and the substitution
is recorded rather than silent. `.oracle_degenerate()` no longer excuses the cell, so it is a real
test again.
**Lives in.** `R/pipeline.R::predict_test`, `R/diagnostics.R::.oracle_degenerate`.

---

## §6. Running it: resources, concurrency and configuration

Not statistics. These are the traps of operating the thing — memory, parallel workers, backups
and config files — and every one of them produced a run that worked while quietly doing the
wrong thing.

### 23. A per-project budget is not a process budget

**Trap.** `dosage_budget_bytes` reads like a memory cap, and the obvious way to use a big-memory
server was to raise it — the README's own example suggested `16e9`. But it caps ONE project's
dense dosage matrix at parse time. `choose_geno_sources` loads *every* project covering a trial at
once, so the process peak is the SUM. On the T3 wheat archive that is about 6x the per-project
figure: two panels alone (project 8160 at 683 x 7.47M, 8124 at 811 x 3.04M) are 20.4 GB and 9.9 GB
dense, and 66 cached projects total 27.9 GB at a 16e9 budget against 11.7 GB at 2e9.
**Symptom.** Latent, because nothing measured it — there was no `gc()`, `object.size` or RSS call
anywhere in the package, and no memory column in `evals`. It surfaced only as the question "can I
run more than one worker?", whose honest answer at `16e9` was no: one worker could hold ~28 GB of
integer dosage, and `.qc_markers` then doubled it by assigning a `double` (`mu` / `round(mu)`) into
an integer matrix, silently coercing the whole thing to 8 bytes per cell on the first imputed
column. Eight of those does not fit in 512 GB.
**Now.** Three changes. `dosage_total_budget_bytes` bounds the sum, serving the largest panels
coarser out of the existing cache (`get_project_dosage`'s `marker_thin`, so no re-download);
`.qc_markers` keeps the matrix integer on the `mean_round` path and subsets once instead of twice,
halving the dominant allocation with bit-identical GRMs; and every evaluation records `peak_r_mb`
so the budget is set from measurement (`tools/report_memory.R`) rather than from arithmetic.
**Lives in.** `R/pipeline.R::.dosage_thin_plan` / `.qc_markers`, `R/memory.R`, `tools/report_memory.R`,
`tools/watch_memory.sh`.

### 24. "SQLite allows one writer" is true and was still the wrong conclusion

**Trap.** The README told you not to point parallel jobs at one `evals.sqlite` — lock contention,
corruption — and the code agreed: `open_store` was a bare `dbConnect` with no WAL and a **zero**
busy timeout, so a second writer got `SQLITE_BUSY` immediately rather than waiting.
**Symptom.** Throughput was capped at one evaluation at a time on a machine with cores to spare,
for a workload that writes one small row per multi-minute evaluation — contention that barely
exists. The real obstacles were elsewhere and unstated: an unguarded `ALTER TABLE` migration loop
that kills whichever worker loses the startup race, `.cache_save` writing `saveRDS` straight to the
final path (so a reader, or the `dosage_*` glob, can hit a half-written file), and
`.ensure_project_vcf` unlinking every copy before downloading — two workers on one project delete
each other's in-flight multi-GB VCF and neither finishes.
**Now.** WAL + a 60 s busy timeout, migrations that tolerate losing the race, atomic
cache writes (temp + `rename`), and a per-project download lock where the loser waits for the
winner's result instead of repeating the download. The remaining constraint is real and worth
stating plainly: **WAL cannot work over NFS**, so the live store must be on local disk, backed up
via `VACUUM INTO` (by any worker -- see #25).
**Lives in.** `R/store.R`, `R/data_access.R::.with_cache_lock` / `.cache_save`, `run_workers.sh`,
`tests/test_concurrency.R`.

### 25. A periodic task at the bottom of a long loop has the loop's period, not its own

**Trap.** `db_backup_minutes = 30` reads like "back up every 30 minutes". The call sat at the
bottom of the evaluation loop and was gated on `leader`, so its real period was
`max(30 min, duration of worker 1's current evaluation)`. With evaluations running 30 min to 33 h,
those are different numbers by three orders of magnitude.
**Symptom.** Found by hand on 2026-08-05: the `/workdir` store's backup was **24 h stale** while
worker 1 sat inside one `em_combine` evaluation. Every completed evaluation since existed only in
the WAL, on the disk the backup exists to protect against losing. Nothing logged a problem,
because from the code's point of view nothing had gone wrong -- the condition simply had not been
reached. The same trap had already been found and fixed for `report.md` (the comment above
`write_report` in `run_optimizer.R` records it) and the reasoning was not carried across to the
backup, which is the file whose loss actually costs work.
**Now.** There is no interval: the backup runs after **every** evaluation, from whichever worker
just finished, and never consults `is_leader`. That is the honest reading of the trap -- an
interval at the bottom of this loop was only ever a floor, and `VACUUM INTO` was measured at
**0.01 s** on stores of 500 to 5,000 rows, so the floor was buying nothing while leaving up to
one un-backed evaluation per worker whenever they all entered long evaluations together. The
final on-exit backup is likewise every-worker: a leader-only one does nothing if worker 1 is the
process that gets OOM-killed. The report states how many rows the backup is **behind** rather
than an age against an interval, because `write_report` is *not* leader-gated and so keeps
updating while a long evaluation runs -- it is the one artefact that can report the stall while
the stall is happening.

The cache sync kept its interval. An rsync tree-walk over thousands of files and gigabytes is
not in the same cost class as a VACUUM, so there the floor earns its place.
**Generalise.** Any periodic work placed at the bottom of a loop inherits the loop's period. If
the loop body can run longer than the interval, the interval is a floor and nothing more; and if
the work is also gated on one process, one slow process disables it entirely. Ask what the *worst*
iteration costs, not the median.

**And it recurred, which is the part worth remembering.** The fix above was applied to the store
backup only. The **cache** backup had the identical shape — `leader &&` at the bottom of the same
loop, plus an `on.exit` flush registered *inside* the `if (leader)` block — and was left in place
for another four days, until the user asked why the reasoning did not apply to it equally. When a
defect is a *shape* rather than a line, grep for the shape: every periodic, leader-gated call in
the same loop was suspect the moment the first one was.

The cache needed one thing the store did not. A `VACUUM INTO` of a 250 KB file is harmless N
times over, so `should_backup_now()` needs no mutual exclusion; an rsync over thousands of files
is not, and dropping the leader gate alone would have traded a stall for a stampede. So
`sync_cache_to_backup()` throttles on a **stamp file** in `cache_backup_dir` and touches it
*before* transferring — claiming afterwards leaves the whole transfer as a window in which every
other worker also starts one. `tests/test_concurrency.R` §6 pins that ordering: with the claim
moved after the rsync, three concurrent workers produce three transfers instead of one.

**What the leader still does, and why it is different.** Restoring the cache at startup stays
leader-only. That is not a throughput throttle but a correctness one: the others must not read a
half-populated tree, which is what `cache_ready_file` exists to signal.

**Then the fix removed the need for a second mechanism.** Once both backups were every-worker
and on a 30-minute interval, `max_hours` — set below the SLURM `--time` so the loop could stop
"cleanly" before the scheduler killed it — stopped earning its place, and was dropped from the
cluster config. A hard kill now costs at most one backup interval and corrupts nothing
(tmp-then-rename backups; WAL is crash-safe). Against that, a margin idles every worker for its
duration and still cannot save the in-flight evaluation, which is only ever checked *between*
evaluations. It had also been silently wrong in the one place it was meant to be exercised:
`--qos=debug --time=00:30:00` left `max_hours = 480`, so the shakeout never rehearsed the
shutdown at all. **A safety net that duplicates a stronger one is not free — it is another
number to keep in step, and it was already out of step.**

*Considered and rejected:* `#SBATCH --signal=B:USR1@300` fires five minutes before the wall
clock and could flush with no idle margin. It needs signal handling in a non-interactive
`Rscript` under `bash … wait`, and it buys back the same ≤30 minutes the periodic backups
already bound. Not worth the machinery.
**Lives in.** `R/store.R::should_backup_now` / `backup_age_minutes`,
`R/data_access.R::sync_cache_to_backup`, `run_optimizer.R`, `R/report.R`,
`tests/test_concurrency.R` §5-6.

### 26. A gitignored config file is frozen at the moment you copied it

**Trap.** The `.local` pattern -- tracked `X.example`, gitignored `X` -- keeps `git pull` from
conflicting on machine-specific values. It is right for that. But `git pull` also never
*updates* the copy, so anything mechanical living in it is pinned to the day it was made.
**Symptom.** After the log directory moved to `$OPTIMIZER_HOME/logs`, a shakeout wrote
`run_wXX.out` to the new place and `slurm-<jid>.out` to the old one. `run_workers.sh` is
tracked, so its half of the change arrived; `container/submit.local.sh` is not, and the copy in
use predated the `--output`/`--chdir` flags, so SLURM fell back to a relative `#SBATCH`
directive. No conflict, no warning, a run that worked -- just an older one. In two days that
one file had drifted four times over.
**Now.** The templates hold values only and source their mechanism from tracked files:
`submit.local.sh` (account, sizes) calls `submit_optimizer()` from `container/lib_submit.sh`;
`settings.local.R` (`dosage_budget_bytes`) splices in `cluster_scratch_paths()` from
`settings.R`. `cluster_scratch_paths()` has to sit in `settings.R` rather than `R/`, because
`.local_overrides()` reads `settings.local.R` before anything in `R/` is sourced. As a
backstop, `lib_submit.sh` exports `OPTIMIZER_SUBMIT_LIB=1` and `t3opt_ceres.sh` warns when it
is absent -- a stale copy still runs, it just says so.
**Generalise.** Ask of every gitignored file: *if I fix a bug in the thing this was copied
from, does the fix reach the copy?* If not, only values belong in it. The tell is a template
that keeps needing edits -- that is not a template, it is code living in the wrong place.
**Lives in.** `container/lib_submit.sh`, `container/submit.local.sh.example`,
`settings.R::cluster_scratch_paths`, `container/settings.local.R.scinet`,
`container/t3opt_ceres.sh`.

### 28. Evaluations differ by orders of magnitude, and nothing bounds them

**Trap.** Reasoning about the run from a typical evaluation -- planning a wall clock, a worker
count or a node size around "about half an hour each".
**The spread, over 121 real evaluations.** Median **29 min**, mean **122 min**, max **33 hours**
-- a 68x range, and the mean is four times the median because the tail dominates. Memory
likewise: `peak_r_mb` median **19 GB**, max **82 GB**, and that is R's heap peak, an
*under*-estimate of true RSS. The whole set cost 246.8 hours of compute. Kernel medians differ
much less than the tail does (`em_combine` 36.4 min, `rkhs_gaussian` 36.1, `vanRaden_single`
22.4), so the method alone does not tell you which evaluations will be the expensive ones --
trial size and panel coverage do.
**Nothing constrains this, deliberately.** There was a `max_eval_minutes` cap; it was removed
because a wall-clock cap **censors non-randomly**. The slow configurations are the thorough
ones, which may also be the best; a capped evaluation stores no score, so a configuration that
always exceeds the cap can never become incumbent however good it is. Set one and the search
"discovers" that cheap pipelines win -- an artifact. Memory is bounded instead, and only where
bounding is safe: `dosage_total_budget_bytes` thins at serve time rather than discarding work.
**So watch it rather than cap it.** `tools/report_timing.R` and `tools/report_memory.R` exist for this.
Size a machine from the *maximum* per worker, not the median -- SLURM kills a job that exceeds
its allocation rather than throttling it -- and expect a job's final hours to lose whatever is
in flight, because a 33-hour evaluation cannot be finished by any amount of margin.
**Lives in.** `tools/report_timing.R`, `tools/report_memory.R`, `settings$dosage_total_budget_bytes`.

### 29. A stop file on durable storage outlives the job that consumed it

**Trap.** `touch "$STOP_FILE"` halts a run cleanly, so the same file was left to be cleaned up
"at startup" by the leader -- `unlink(settings$stop_file)` in `run_optimizer.R`. But the state
directory is durable by design, so the file survives the job, and the *launcher*
`run_workers.sh` had its own pre-flight check that refused to start while it existed. The check
ran in a different process, in a different language, 38 lines before the first worker was
spawned, so the clearing code was unreachable on every path that went through the launcher --
which is every SLURM launch and every multi-worker interactive one. Stopping a job through the
stop file therefore bricked the next one: `exit 1`, no evaluations, and `--dependency=singleton`
released the job after that to fail identically. The runbook confidently documented the
behaviour that could not happen.
**Why the tests did not catch it.** The oracle (`tests/test_subtasks.R`, "a stale stop-file does
not halt the next run") calls `run_optimizer()` directly in R. It asserts exactly the right
thing about exactly the path nobody uses in production.
**Generalise.** When a guard and its remedy live in different processes, work out which one
runs first -- an unreachable remedy reads like working code. And a cleanup step belongs in
whatever *starts* the thing, not in the thing being started: only the launcher runs before all
N workers, so only the launcher can clear a flag that all N of them test.
**Lives in.** `run_workers.sh` (clear-or-refuse, keyed on `OPTIMIZER_FIRST_WORKER`),
`run_optimizer.R` (still clears it, for a bare `Rscript run_optimizer.R`).

### 30. Catching `warning` alongside `error` threw away a good fit, and took contender replication with it

**Trap.** `aggregate_scores` guarded the random-effects fit with
`tryCatch(.blup_scores(evals), error = function(e) NULL, warning = function(w) NULL)`. The
warning handler reads as prudence and is not: `lmer` returns a perfectly usable fit and *then*
warns, and this design warns routinely, because `sd_config` (~0.04) sits near the boundary
against `sd_trial` (~0.10) — the regime that produces "Model failed to converge". Every such
fit was discarded and silently replaced by the pooled estimate.
**What it cost, which was not the obvious thing.** The report losing `sd_trial`/`sd_config`/
`sd_resid` was the visible symptom. The real damage: `se` is supplied *only* by that fit, and
`.contenders()` filters on `is.finite(se)`, so no fit meant no contenders, which meant the
`$extra` tier of the replication backlog was permanently empty and `replicate_every` had
nothing to ration. A whole scheduling tier was switched off by an error handler. Measured on a
4-worker simulate run: `estimator: pooled, contenders: 0` before, `blup, contenders: 8` after.
**Generalise.** `tryCatch(warning=)` does not "handle" a warning, it *aborts the call* — the
computation never returns. Use `withCallingHandlers` + `invokeRestart("muffleWarning")` to keep
the result and record the message. And never let a fallback be silent: every path that
substitutes a weaker estimate now writes `estimator_note`, which the report prints, because
"pooled" with nothing after it looked like a choice rather than a failure.
**Lives in.** `R/optimizer.R::aggregate_scores` / `.blup_skip_reason`, `R/report.R`.

### 31. Excluding in-flight work needs the exclusion where the write is, not where the decision is

**Trap.** Two workers could run the identical `(config, trial, scheme)` — deterministic, so the
second recomputes a known number and double-counts it in that config's mean. `choose_trial`
excluded trials from `read_evals()`, which sees only *committed* rows, so a one-entry backlog
handed every worker the same trial. Adding a `claims` table with `PRIMARY KEY (config_hash,
trial_id, scheme)` fixed the concurrent case cleanly — and duplicates kept appearing.
**The half-fix.** The claim answered "is anyone running this?" but not "is this already done?".
The holder commits its row and releases ~2 ms later; a worker whose selection read the store
just before that commit finds the pair neither claimed nor (in its stale copy) evaluated, and
claims it legitimately. Six duplicates in 240 evaluations, all pairs of different workers in
the same second. Widening the advisory filter cannot close this: **any** read is stale by the
time the caller acts on it.
**Generalise.** A test-and-set must test everything that makes the action wrong, in the same
statement that performs it. Here that is one `INSERT ... SELECT ... WHERE NOT EXISTS (SELECT 1
FROM evals ...)`: in flight *or* already finished, decided at the write. The advisory filter
still earns its place — it stops most collisions before an expensive selection — but it is an
optimization, never the guarantee.
**And staleness is a dead owner, never elapsed time.** A claim held for three days is a
three-day evaluation (33 h observed maximum, #28); expiring it would start a duplicate of live
work. `.pid_alive()` is the test, run only when a claim is refused, so the `ps` sweep is paid
for on contention rather than on every selection.
**Lives in.** `R/store.R::claim_eval` / `.purge_dead_claims` / `.pid_alive`,
`run_optimizer.R::optimizer_step`, `R/optimizer.R::choose_trial`.

### 32. A backup that is a copy of the working store keeps only what the restore brought back

**Trap.** A job on a fresh node starts with an empty node-local store, and
`restore_store_from_backup()` copied `evals` alone. The next `backup_store()` is a `VACUUM INTO`
of that store and replaces the backup file. So every restart on a new node dropped the `runs`
records of all past domains, while their evaluations survived. Nothing noticed until
`dev/holdout_test.R` was asked for a past domain: its run id was correct, and the store holding
38,340 evaluations had only the current run's record.
**Side effect.** Losing the record also meant each chained job pinned its trial universe again
from the catalogue, so the protection against a mid-run trial rename (the August stall) lasted
one job, not one run.
**Generalise.** When a backup is a snapshot of the working copy rather than an accumulating
archive, the restore must bring back *every* table worth keeping. Anything it skips is deleted
at the next backup.
**Recovery.** A record already lost cannot be rebuilt from `evals`; `dev/holdout_test.R`
resolves the domain from today's catalogue when no record matches (`resolve_trial_universe()`).
**Lives in.** `R/store.R::restore_store_from_backup`, `R/data_access.R::resolve_trial_universe`,
`dev/holdout_test.R`.

### 34. One quantity, two thresholds: the band between them was evaluations paid for and thrown away

**Trap.** Gating *feasibility* and *usability* on the same quantity with different numbers.
`score_predictions()` produced a score from as few as **5** joined, finite accession pairs — the
Predictathon's own "task can fail" rule — while `aggregate_scores()` admitted a score to the
estimates only at **10** or more (`min_n_test`, a hardcoded default argument, not a setting).
Each number was defensible alone. Together they defined a band in which an evaluation ran to
completion, cost a full pipeline, was stored with `status = "ok"` and a real correlation, and
was then invisible to the BLUP, the contender rule, the elite pool and `n_ok`.
**Symptom.** Nothing. Both gates behaved exactly as written; no warning, no failure row, and
the report counted those evaluations as successes. The waste was only visible by crossing
`n_test` with `trial_id`: on run `9a2664962c01`, six of the fourteen trials ever evaluated
returned `n_test = 7` on **every** evaluation, so **2029 of 5631 successful evaluations (36%)**
were discarded by every estimate, and the maximum attainable replication for any configuration
in that run was 8 usable evaluations rather than 15. The same crossing on `80cfbd7c27e9` loses
10% and on `46b0d4b39c79` none, so the cost depends entirely on the target domain — which is
why no single run's health said anything about it.
**Aggravator.** The Fisher-z weighting (`n_test - 3`) made the pair *look* harmless: a short
correlation already enters the pooled mean at near-zero weight, so the hard floor reads like a
belt-and-braces refinement. It is not, because the surrogate is **not** weighted — see #35 —
and because a weight cannot refund an evaluation that should never have been attempted.
**Now.** One setting, `settings$min_test_acc`, default **10**, read through `.min_test()` and
applied at both ends: `score_predictions()` declines the task below it, so no row is created,
and `.fisher_weights()` marks a row `.usable` only at or above it. `min_test` is a required
argument of both `score_predictions()` and `.fisher_weights()`, and `aggregate_scores()` takes
`settings` — there is no local default left to drift. The `.usable` test is still live code
because the archive holds rows written under the old pair.
**Cost of the fix, stated plainly.** 10 declines tasks the Predictathon would score at 5. That
is deliberate: a correlation on 5 accessions is not worth a 29-minute median evaluation on
either side of the gate, and the alternative — scoring it and then ignoring it — is what this
entry is about.
**Generalise.** When two stages of a pipeline test the same quantity, they must read the same
constant from the same place, and that place should be the settings file. Two independently
reasonable thresholds on one quantity do not compose into a reasonable policy; they compose
into a silent dead zone whose width nobody chose and whose cost nobody measured. The smell to
look for is a hardcoded default argument that shadows a real setting.
**Lives in.** `settings.R::min_test_acc`, `R/pipeline.R::.min_test`,
`R/evaluate.R::score_predictions`, `R/optimizer.R::.fisher_weights` / `aggregate_scores`.

### 35. The surrogate sees raw, unweighted correlations — not the z-transformed weighted ones

**Trap.** Assuming that because `aggregate_scores()` does the statistically careful thing —
Fisher-z transform, weight `n_test - 3`, trial effect blocked out by a mixed model — everything
downstream of a score does too. It is a natural assumption and it is what makes #34's two
thresholds look safe.
**Symptom.** None visible. The surrogate predicts, the search proceeds. But `choose_config()`
builds the forest's training rows with `dplyr::filter(evals, is.finite(score))`
(`R/optimizer.R:521`) — `is.finite`, not `.usable` — and passes raw `score` to
`fit_surrogate_merf()`, which applies **no** weights, **no** z-transform, and **no** `n_test`
floor at any point: not in the `lmer` call that estimates the trial effect, not in the `rpart`
trees. So a 5-accession correlation of 0.9 and a 200-accession correlation of 0.9 are the same
observation to the forest, and rows the BLUP refuses to use still shape every prediction and
every EI value.
**Now.** Closed, 2026-11-03. `choose_config()` runs its slice through `.fisher_weights()` — the
estimator's own helper — keeps only `.usable` rows, fits on `.z`, and passes `.w` as `rpart`
case weights and as `lmer` weights inside `fit_surrogate_merf()`. One floor, one weight, one
transform, one helper, two consumers. The cap on that weight is #36; putting `best` on the same
scale as the predictions took a second step, and is #37.
**Generalise.** Two consumers of the same column are two chances to disagree about what the
column means. Being careful in one estimator does not propagate; the care has to be applied
where each consumer reads, or the column has to arrive pre-weighted. The structural fix is not
"weight the second consumer too" but "make both call the same function", which is why
`.fisher_weights()` is now the only place a row's worth is decided.
**Lives in.** `R/optimizer.R::choose_config` (the row filter),
`R/surrogate.R::fit_surrogate_merf` / `fit_surrogate`.

### 36. `n_test - 3` is the wrong weight, because sampling error is not the only error

**Trap.** Taking the textbook inverse-variance weight for a correlation, `n_test - 3`, as
correct because the Fisher transform makes `Var(z) = 1/(n-3)` exactly. It is the right weight
for *resampling the test set*. This design never resamples a test set: every configuration on a
trial is scored against the **same** accessions, and a configuration can also genuinely perform
differently on that trial than on an average one. That second error — config × trial interaction
— does **not** shrink as the trial grows. A bigger trial measures the trial's idiosyncrasy more
precisely; it does not make the idiosyncrasy smaller.
**Symptom.** Weights spanning 7 to 1434 within one run (205×), so a handful of large trials
decided every split and every pooled mean. And a calibration failure hiding in plain sight:
`lmer` reports `sd_resid` at unit weight, so if `n_test - 3` were the whole variance story the
fitted `sd_resid` would be **1**. Measured across seven runs it spans **0.096 to 4.72** — a
factor of 49, in *both* directions. Above 1 is the unmodelled interaction; below 1 is the shared
test set, whose realised sampling error is common to every configuration on that trial and is
therefore absorbed by `(1 | trial_id)` rather than appearing in the residual.
**Now.** Two precisions in series, `1/w = 1/(n_test - 3) + 1/weight_cap`, in `.eval_weight()`,
used by `.fisher_weights()` and so by both the aggregate and the surrogate. `weight_cap` is the
weight an unboundedly large trial gets, so it reads as a trial size: the default **50** means
"no trial counts for more than one with 53 accessions", and collapses the within-run weight span
from 205× to 8×.
**What the calibration actually said, and why the default is not it.** Solving `sd_resid = 1`
per run gives `C ≈ 33` and `41` for the two CV0 runs and `C ≈ 608` and `1058` for the two CV00
runs — tight within scheme, 25× apart across it, which fits the mechanism (CV00 masks the focal
genotypes, so configurations are forced onto more similar information and interact less with the
trial). The other three runs cannot inform it at all: their `n_test` is near-constant (127–130),
which makes the fit *invariant* to the cap — identical `sd_trial`, `sd_config` and `n*` at every
grid point. A single conservative 50 was chosen over two per-scheme numbers fitted on two runs
each. It is right for CV0 and ~20× too tight for CV00.
**Cost of that choice, stated plainly.** Capping at 50 drives `sd_config` to exactly **0** in
both 24-trial runs (`80cfbd7c27e9`, `e7517df121d5`), where uncapped they fitted 0.0135 and
0.0305. `80cfbd7c27e9` keeps a non-zero 0.0105 at its own calibrated 1058, so for that run the
collapse is the cap's doing, not the data's. And `sd_config = 0` is worse than it sounds: every
BLUP shrinks to exactly the grand mean, so all configurations **tie** at the maximum, all of them
pass the contender test (`in_play` goes to 2147/2147 and 705/705), and the top-8 cap hands
`get_elites` an arbitrary eight. The ranking stops carrying information rather than merely
getting noisier.
**But that zero is singleton dilution, not a lost signal — measured.** Re-fitting
`80cfbd7c27e9` on replicated configurations only takes `sd_config` from 0.0000 to **0.0068**
(`n_ok >= 4`) and **0.0134** (`n_ok >= 6`) — back to its uncapped whole-population value. The
capped fit is not blind; it is telling the truth about a population that is mostly singletons,
and `config_replication` plus contender replication are the mechanisms that answer it. The
contrast is the evidence: on `e7517df121d5` and `f92e28ec6415`, `sd_config` is 0 at **every**
replication threshold, so there the zero is the answer. Which is why one conservative cap was
acceptable after all — on a 24-trial CV00 target, the shape of the live run, it costs
replication effort, not information.
**Generalise.** A variance model inherited from a textbook is a hypothesis about your design,
not a property of your statistic. It is testable almost for free whenever the fit reports a
residual at unit weight: calibrated weights put that residual at 1, and any other value is the
size of what the weight does not model. Check it before trusting the weight, and check it per
run — ours disagreed by 49×.
**Lives in.** `settings.R::weight_cap`, `R/optimizer.R::.eval_weight` / `.fisher_weights`,
`R/surrogate.R::fit_surrogate` / `fit_surrogate_merf`, `dev/calibrate_weights.R`.

### 37. EI compared a shrunken incumbent against unshrunken predictions, and that killed exploration

**Trap.** Reaching for the best number you have as the EI reference point. `max(mean_score)` is
the best estimate of the best configuration's quality, so it looks like the right incumbent —
and it is, as an estimate. It is the wrong *reference*, because EI does not compare estimates,
it compares a prediction against a reference that must be produced the same way.
**Symptom.** Acquisition quietly stops exploring. Nothing errors, nothing warns, EI is still
computed and still has a well-defined argmax. But `best <- max(.to_z(scored$mean_score))` is a
BLUP, shrunk toward the grand mean in proportion to `sd_config` (SEARCH_ALGORITHM §4.5), while
`pr$mean` is a forest prediction with no shrinkage at all. Write EI out:

```
EI = (m - best - xi)·Phi(z) + sd·phi(z),      z = (m - best - xi)/sd
```

Driving `best` down sends `z -> +inf`, so `Phi(z) -> 1`, `phi(z) -> 0`, and `EI -> m - best - xi`.
The ranking does not flatten — it **collapses onto the predicted mean**. `sd·phi(z)` is the only
term that ever lets EI prefer a candidate the forest is unsure about, so losing it converts
Expected Improvement into plain greedy exploitation of the surrogate mean. And the error grows
with exactly the thing that should trigger more caution: the less the configurations are
distinguishable, the smaller `sd_config`, the harder the shrinkage, the greedier the search. In
the limiting case `sd_config` fits to 0, every BLUP *is* the grand mean, and `best` was the
grand mean of the whole run.
**Now.** Closed, 2026-11-03. `best` is the forest's own trial-marginal prediction at the
features of the configurations already evaluated, through the same `.marginal()` closure that
predicts the candidates — one estimator, one scale, by construction. **Measured** on the final
rows of four runs at `weight_cap = 50`: `best` rises 0.1211 -> 0.2239 (`80cfbd7c27e9`) and
0.3516 -> 0.5228 (`e7517df121d5`), the two runs where `sd_config` fits to 0, and Spearman(EI,
predicted mean) falls from 0.682 to 0.260 and from 0.766 to -0.017. Where shrinkage is mild
(`34e15bd03838`, `sd_config = 0.1001`) `best` moves -0.02 and little changes. The size of the
correction tracks the shrinkage, which is the algebra's own prediction.
**Why the repair only runs one way.** The symmetric-looking alternative — shrink `m` to the
posterior scale instead of deregressing `best` — is not available. A candidate's shrinkage
factor is `lambda_i = W_i/sd_resid² / (1/sd_config² + W_i/sd_resid²)` and an unevaluated
candidate has `W_i = 0`, so `lambda = 0` and every candidate would collapse to `mu_hat`. The
posterior scale has no room for a configuration with no data. Harmonisation has to go toward the
unshrunken scale, which is why "use the model's own prediction at the evaluated points" is the
standard EI reference under noisy observations rather than a local trick.
**What this does NOT fix, measured.** On `9a2664962c01` the change does nothing at all —
Spearman(EI, mean) stays at exactly 1.000 — because the uncertainty floor
`sd_floor <- 0.25 * sd(y_obs)` binds for **129 of 129** candidates, and under a uniform `sd` EI
is monotone in `m` whatever `best` is. There are two independent routes to the same degeneracy
and only one is closed. The floor is the open one.
**Open: the reference is partly in-sample.** Predicting at a configuration's own features uses
trees that saw its rows. Measured against the unshrunk pooled means the forest smooths far
harder than it memorises — `max(prediction)` is 0.22-0.52 where `max(pooled z)` is 0.70-0.93,
and restricting to `n_ok >= 3` barely moves it — and the bias errs toward a higher `best`, hence
more weight on `sd`, which is the safe direction here. An out-of-bag reference would remove it
at the cost of storing per-tree in-bag indices; not built.
**Generalise.** Whenever two numbers are subtracted, ask what estimator produced each. Shrinkage
is the usual culprit because it is invisible at the call site: both sides are "the score",
on the same units, back-transformed the same way, and only the *estimation* scale differs. The
defence is structural — make both sides come out of the same function, as `.marginal()` now
does — not a correction factor, which would need the very variance component that was unstable.
**Lives in.** `R/optimizer.R::choose_config` (the `.marginal()` closure and `best`),
`R/surrogate.R::expected_improvement`, `docs/SEARCH_ALGORITHM.md` §8.3-8.4.

### 38. The log knew which phase proposed a configuration; the store did not

**Trap.** Logging a fact and recording it are not the same thing, and the difference is invisible
until someone asks a question about the past. `choose_config()` has always returned `source` --
`seed:Prediction3`, `random_init`, `acquisition_merf`, `random_fallback` -- and the main loop has
always printed it (`src=acquisition` in every iteration line). It was never written to `evals`.
**Symptom.** A question that should have been a one-line query turned out to be unanswerable.
`sd_config` is estimated over the configurations *in the fit*, and the search narrows that
population deliberately, so the estimate falls as the search concentrates (#33, SEARCH_ALGORITHM
§11). The fix is to estimate it instead over the configurations chosen **independently of the
surrogate** -- the seeds and the Phase 3 random draws -- which is an unbiased sample of the
universe and whose rows are sitting in the store already. Except that nothing in `evals` says
which rows those are. The best available reconstruction is "the first `n_random_init + 5`
distinct `config_hash` by `ts`", which is approximate (concurrent workers interleave), fragile,
and silently wrong for any run whose phases did not proceed in the expected order.
**Measured, on that reconstruction.** Where the reconstructed set is replicated enough to
identify a variance component at all, it says what the theory predicts: on `80cfbd7c27e9`
(median `n_ok` 5) `sd_config` is **0.0130** over the selection-independent set against 0.0000
over the whole population, and 0.0130 agrees with both the uncapped whole-population fit
(0.0135) and the `n_ok >= 6` restriction (0.0134) -- three routes to the same number. On
`9a2664962c01` (median 6) it is 0.0246 against 0.0078. On the three runs whose reconstructed set
has median `n_ok` of 1-3 the estimate is noise and lands on both sides. So the idea works and the
archive can only half-demonstrate it.
**Now.** `evals` has a `source TEXT` column, added by the same `ALTER TABLE` migration loop that
added `build` and `run_id`, so old stores gain it on open. `store_eval()` takes `source` and
`optimizer_step()` passes `choice$source`. Rows written before 2026-11-03 are `NA` and stay
unclassifiable -- the column is not retroactive, which is the whole point of this entry.
**Generalise.** When a decision procedure has branches, record which branch fired, in the same
row as the outcome. It costs one column and it is the only way to ask afterwards whether a
subset of the data was selected on the thing you are now trying to estimate. The tell is a
variable that exists in the code, is interesting enough to print, and is not in the schema: the
log is a stream nobody can join against, and by the time the question arrives the stream is
gone. Selection provenance in particular is not recoverable by inference -- `ts` ordering is a
guess, not a record.
**Lives in.** `R/store.R` (schema, the `add` migration vector, `store_eval`),
`run_optimizer.R::optimizer_step`.

### 39. Two unrelated quantities shared one setting, and the coupling was load-bearing

**Trap.** `min_obs` -- the rows a forest needs before it is worth fitting -- was passed as
`settings$n_random_init` at all three `fit_surrogate*` call sites. The two numbers want to be
equal at the start of a run, so one setting looked like economy rather than conflation. They are
different questions: `n_random_init` is a search-policy budget ("explore at random until this
many configurations are scored"), `min_obs` is a statistical floor on a regression fit.
**Symptom.** Neither could be changed alone. Raising the random-exploration budget to 100 also
demanded 100 rows before any surrogate would fit; lowering the fit threshold also cut exploration
short. Worse, `n_random_init` was at the time one of the ten run-defining keys and `min_obs` has
no business being one, so a purely statistical choice could not be made without rotating the
run's identity. (Both halves of that have since moved: `run_identity = "target"` makes
`n_random_init` policy rather than identifying, and #45 records that rotating `run_id` never hid
the rows in the first place. SEARCH_ALGORITHM §13.)
**What the coupling was silently doing.** Decoupling them broke a test, which is how the load-
bearing part surfaced: `min_obs` counts the rows of whichever mode is being tried, and the modes
count differently -- MERF and BLOCKED get one row per `(config, trial)`, POOLED gets one row per
*configuration*. Phase 3 guarantees exactly `n_random_init` scored configurations before
acquisition begins, so `min_obs = n_random_init` was precisely the largest threshold at which
POOLED could ever fit. Setting the new `surrogate_min_obs` above `n_random_init` therefore
retires the POOLED fallback without saying so; the search falls through to `random_fallback`
instead, which is safe but spends the step.
**Now.** `settings$surrogate_min_obs` (default 25, equal to `n_random_init`, so behaviour is
unchanged), read through `.surrogate_min_obs()`, validated as an integer >= 1, and not
run-defining. The mode-dependent row count and the POOLED consequence are written next to the
setting and in SEARCH_ALGORITHM §8.2.
**Generalise.** Sharing one setting between two quantities that merely happen to be equal is a
latent constraint, and removing it can remove a guarantee nobody wrote down. Before splitting a
coupled constant, ask what the equality was protecting -- the answer is often in a test that
breaks for a reason unrelated to the change. The converse trap is #34: one quantity read from
two places. Both come from the same failure to ask "is this the same number, or two numbers that
agree today?"
**Lives in.** `settings.R::surrogate_min_obs`, `R/optimizer.R::.surrogate_min_obs` and
`choose_config`, `dev/measure_churn.R`.

### 40. `sd_config` was measured on the population the search had already narrowed

**Symptom.** The shrinkage prior tightened *because the search was working.* `sd_config` is the SD
of the configuration effects in whatever population is currently in the fit, and by late in a run
most evaluated configurations are recombinations of a few parents. On `9a2664962c01` the fitted
value fell from 0.0251 to 0.0078 over 6236 evaluations; on `80cfbd7c27e9` and `f92e28ec6415` it
reached the boundary at 0, where every BLUP ties at the grand mean, every `se` is 0, and the
contender rule is inoperative at every `z`. Nothing errored, and the number was a correct REML
estimate of the quantity it was actually estimating.

**Why it mattered three times over.** That one number is the shrinkage prior, it is the unit the
contender test is denominated in (`se <= sd_config`), and it feeds back into the homogenization
that lowered it. So a search that concentrated became a search that could no longer tell anything
apart, and the mechanism looked like success.

**What made the repair possible.** #38 — storing `evals.source`. Once the store records which
phase of `choose_config()` proposed a configuration, the rows drawn *without* reference to the
surrogate can be named, and `sd_config` can be estimated on those alone. That set is a sample of
the configuration space; the full evaluated set is a sample of the search's trajectory. Seeds are
excluded: five hand-picked submissions are not a sample of anything.

**Measured.** Replaying the real runs (`dev/sweep_contender_rule.R`), the reference estimate is
flat across a run where the fitted one collapses — `9a2664962c01`, fitted 0.0251 → 0.0078,
reference 0.0240 → 0.0246 — and it exists where the fitted one is zero: `80cfbd7c27e9` 0 → 0.0130,
`e7517df121d5` 0 → 0.0517, `f92e28ec6415` 0 → 0.0250. On `34e15bd03838` it runs the *other* way
(fitted 0.100, reference 0.056), which is the honest counter-case: a run whose evaluated set
genuinely is more spread out than the space gets a *tighter* prior, not a looser one. The
consequence for the search is large: median in-play share goes from 100% to 0.7% on
`80cfbd7c27e9` and from 100% to 0.05% on `f92e28ec6415` — from a set containing every scored
configuration to one the cap of 8 can actually serve.

**Two things the repair does not fix.** The reference population is frozen when Phase 3 ends, so
its estimate stops improving — with `n_random_init = 25` and `config_replication = 2` it rests on
about 2 rows per configuration, and on `7d75f214584a` that is thin enough that the reference fit
itself returns a boundary 0 and the code declines (correctly: a pin of 0 is strictly worse than the
fitted value it would replace). A persistent exploration rate, or a larger `n_random_init`, is the
lever, and neither is built. Second, the pinned path re-shrinks analytically while the unpinned
path uses lme4's `postVar`, so switching `sd_config_source` switches the `se` *estimator* as well
as the prior; the analytic form runs a few percent low, which is immaterial beside the 2–4x prior
differences above but is a real confound in any small comparison.

**Generalise.** When an estimate is used to *decide* what to sample next, ask what population it
is an estimate *of*. A quantity computed from the search's own output describes the search, and
feeding it back as a prior closes a loop that looks like convergence from the inside. The fix is
not a better estimator — it is a sample that the decision did not choose.

**Lives in.** `R/optimizer.R::.reference_hashes`, `.reference_sd_config`, `.blup_scores`
(`sd_config_pin`), `settings.R::sd_config_source`/`reference_min_configs`,
`dev/sweep_contender_rule.R`, SEARCH_ALGORITHM §11.

### 41. The contender bar was a point estimate, so a lucky singleton could evict the leader

**Symptom.** The in-play set did not drift over a run, it *thrashed*: on `9a2664962c01` at
`contender_z = 0.5` the count went 5, 33, 4, 66, 7, 24, 24, 156, 120, 131, 144, 99, 156, 231, 254,
88, 90, 259, 271, 84, 61, 156, 7, 7, 13 across 25 checkpoints while `sd_config` declined smoothly.
Membership in a set that governs *all* further replication (§7: dropping out costs everything) was
turning over almost completely between refits, and there was no stage signal for an adaptive `z` to
track.

**Cause.** The rule compared an optimistic bound against a *point* estimate,
`mean_i + z*se_i >= max_j(mean_j)`. The bar was therefore the maximum of a set of noisy
estimates — the most upward-biased order statistic in the set — so one thinly replicated
configuration with a lucky draw set the bar for everyone, including for a well-replicated leader it
could evict. The asymmetry is also unprincipled on its own terms: the candidate is charged for its
uncertainty and the leader gets a free pass on its.

**Fix.** The racing comparison, `mean_i + z*se_i >= max_j(mean_j - z*se_j)` — optimistic against
pessimistic, both sides charged the same `z`, which is the criterion a racing algorithm eliminates
on (intervals must stop overlapping). A lucky singleton can now *join* the set but cannot
*redefine* it, because setting the racing bar takes a high mean **and** a small `se`.

**What the measurement actually said, including the part that deflates it.** At a fixed `z` racing
is more stable at every run tested, by 2.4–4.6x at `z = 0.5`. But it also produces a set 3–5x
larger, and larger sets are mechanically steadier, so the fair comparison is at matched set size —
and there the advantage mostly evaporates (two runs of three better, one worse). Once the
reference prior of #40 is in place it evaporates further: at comparable set size racing and the
point rule are a wash on stability. The rule is kept because its failure mode is structurally
worse, not because the stability gain survived. **The prior was the lever; the comparison was
not.** Two independent diagnoses had pointed at the bar, and only measurement separated them.

**And the stability metric itself turned out to reward vacuity.** `jump` falls as the set grows,
so the old settings' apparent steadiness was the steadiness of a set containing 66-100% of the
population — it could not move because it already held everything. On the settings that fixed the
set size, `jump` got *worse* on four of seven runs. Any scale-free turnover measure has to be read
beside the quantity it is measuring; here that was stated as a caveat first and then confirmed as
the dominant effect, which is the right order but was close to being forgotten in between.

**Two consequences of adopting it.** `z` does double duty under racing — it inflates the claim and
lowers the bar — so the set grows about twice as fast in `z`, and the default fell from 1 to 0.25
to hold the operating point. And `mean_score` is on the correlation scale while `se` is on the z
scale, so both bars mix the two; every `z` here is calibrated against that mixed comparison, so
correcting the scale would silently move the operating point and needs a re-sweep, not an edit.

**Lives in.** `R/optimizer.R::.in_play`/`.contender_spec`, `settings.R::contender_rule`/
`contender_z`, `dev/sweep_contender_z.R` (which measured the point rule and is pinned to it),
`dev/sweep_contender_rule.R`, SEARCH_ALGORITHM §7.

### 42. Replication and parenthood asked the same question two different ways

**Symptom.** Two decisions turn on "is this configuration potentially good": whether it earns
further evaluations, and whether it is a parent for crossover and mutation. The first consulted
uncertainty — the contender set, `mean + z*se` against a bar. The second did not: parents were the
top `settings$n_elites` by point estimate, `get_elites()` sorting on `mean_score` with `se`
unread. So a configuration could be judged good enough to deepen and not good enough to breed
from, or the reverse, with nothing anywhere giving a reason for the asymmetry. §4.4's claim that
`se` enters the search in exactly one place was the symptom stated as a feature.

**Fix.** One definition, `.in_play()`, read by both. Parents are the in-play set, capped at
`n_parents` (32) and floored at `n_elites` (8).

**Why both bounds are load-bearing.** The cap: when the estimator cannot separate anything the
in-play set is *every* scored configuration, and unbounded parenthood is then recombination among
all of them — ordered by UCB, so the cap keeps the most promising. The floor: the in-play set is
*empty* whenever the random-effects fit did not run, because no fit means no `se` and `.in_play()`
drops every row without one — and an empty parent pool turns `propose_candidates()` into pure
random draws with nothing saying so. That is #30's failure in a new place.

**A worry that measurement refuted.** UCB ranking was expected to flood the parent pool with
singletons, since `se` is largest where replication is smallest. It does the opposite in the usable
range: median `n_ok` of the UCB top 8 is 4.5–5 at `z <= 0.75` and falls to 2 only at `z >= 1.5`.
Shrinkage is why — a singleton's BLUP is pulled so far toward the grand mean that its inflated
`se` does not recover the gap.

**Generalise.** When two decisions depend on the same latent judgement, look for the one that was
implemented first and ask whether the other inherited its criterion or just never got one. The
asymmetry here was not a design choice anybody defended; it was two functions written at different
times.

**Lives in.** `R/optimizer.R::get_parents`/`get_elites`/`.in_play`,
`settings.R::n_parents`/`n_elites`, SEARCH_ALGORITHM §7 and §8.5.

### 43. The contender bound added a z-scale standard error to a correlation-scale mean

**Symptom.** None visible. `.contenders()` computed `mean_score + z * se` where `mean_score` is
`tanh(mu + c_i)` on the correlation scale and `se` is the posterior SD on the Fisher-z scale, so
the sum was a quantity on neither scale. Nothing errored, every number looked plausible, and at
this project's score range (r ~ 0.15) the two differ by about 2% — which is exactly why it
survived: the report even carried a comment saying the report "does not convert" because the
difference is small.

**Why it still mattered.** Two reasons, and the second is the one that forced the fix.

1. `tanh'(x) = 1 - r²` runs 0.99 at r = 0.1 but 0.79 at r = 0.45 and 0.19 at r = 0.9, so the
   error is not uniform: it *overstates* the optimistic bound of a high-scoring configuration
   relative to a low-scoring one. The bias runs in favour of whatever is already winning.
2. `resolution_delta` is a difference in predictive **accuracy**. That only means something on
   the correlation scale, so adding it to a comparison that was on neither scale had no defined
   interpretation at all. The resolution floor could not be specified honestly until the
   comparison it modifies was on one scale.

**Fix.** `.blup_scores()` now returns `z_mean` beside `blup`, so the estimate and its `se` are
available on the same scale. `.in_play()` forms the interval on z and transforms the
**endpoints**: `ucb_i = tanh(z_mean_i + z*se_i)`, `bar = max_j tanh(z_mean_j - z*se_j)`. Exact,
not a delta-method approximation. `atanh(mean_score)` is the fallback for a caller supplying
only `mean_score`, which is exact for a BLUP row.

**What changed and what did not.** Because `tanh` is monotone, the fix changes **no** in/out
decision for `z` alone — the ordering and the membership are identical. What moves is anything
involving `delta`, and the measured effect on real runs is modest but not nil: median in-play at
the live spec went 25.0 → 23.5, 53.0 → 34.5, 179.0 → 141.5 on three runs. So the calibration of
`contender_z` carried over, and the re-sweep confirmed rather than redid it.

**Generalise.** A transform applied to one side of a comparison and not the other is invisible
when the transform is near-linear over the observed range, and small numerical disagreement is
the easiest thing in the world to write a comment excusing. The tell is not the size of the
error, it is that the two operands have different units. Here the giveaway had been sitting in a
comment in `R/report.R` for months: "`se` here is on the Fisher-z scale ... while mean_score is
back-transformed ... the two differ by ~2%, which is why the report does not convert."

**Lives in.** `R/optimizer.R::.blup_scores` (`z_mean`), `.in_play`, `aggregate_scores`,
SEARCH_ALGORITHM §7.

### 44. A prior estimated from 24 configurations was pinned as if it were known

**Symptom.** The reference-population `sd_config` of #40 fixed the collapse it was built for and
then over-corrected on the thinnest run. `e7517df121d5`'s fitted `sd_config` is 0; its reference
estimate is 0.0517; pinning that resolved its contender set to **2 configurations out of 300**
after 2964 evaluations. The pipeline was working exactly as designed and the answer was absurd.

**Cause.** The reference estimate becomes a *fixed* prior, and nothing anywhere allows for its
own uncertainty. It is fitted on as few as `reference_min_configs` = 20 configurations, and the
floor counts **configurations**, not rows. `e7517df121d5` had 37 usable rows across 24 reference
configurations — 1.5 rows each — which is enough to run a crossed fit and not nearly enough to
estimate a variance component from it.

**How it was caught.** Not by inspection. The per-run table showed two runs resolving to 1–2
configurations in play, which is a strong claim, so the estimate they rested on was resampled:
refit each reference population on two thirds of its configurations, 25 draws. Measured at the
full prefix of each run (`dev/`, the profile behind the table below), the spread of the draws
separates the runs cleanly even where the collapse count does not:

| run | reference `sd_config` | IQR of the draws | Q3/Q1 | collapses | usable rows / configs |
|---|---|---|---|---|---|
| `34e15bd03838` | 0.0560 | 0.0561–0.0657 | 1.2 | 0/25 | 86 / 25 |
| `f92e28ec6415` | 0.0251 | 0.0209–0.0299 | 1.4 | 0/25 | 50 / 25 |
| `9a2664962c01` | 0.0246 | 0.0207–0.0284 | 1.4 | 1/25 | 86 / 24 |
| `46b0d4b39c79` | 0.0273 | 0.0156–0.0333 | 2.1 | 4/25 | 32 / 20 |
| `e7517df121d5` | **0.0517** | **0.0265–0.0970** | **3.7** | 5/25 | **37 / 24** |
| `80cfbd7c27e9` | 0.0130 | 0.0053–0.0188 | 3.5 | 6/25 | 169 / 25 |

`e7517df121d5` has 1.5 usable rows per reference configuration, and its quartiles span a factor
of 3.7. The number is an artefact of which configurations happened to be in it.

**Fix.** `reference_stability_draws` (25) runs that resampling inside `.reference_sd_config()`
on every fit and refuses the estimate when more than a quarter of the draws collapse to the
boundary. A row-count floor would **not** have worked: `80cfbd7c27e9` has 6.8 usable rows per
reference configuration and is the *least* stable of the seven by collapse count, because its 25
configurations are spread over 24 trials and few of them share one. Resampling measures whether
the answer depends on the sample; a row count is a proxy for that, and a bad one here.

**What the gate does and does not catch, measured.** On this store it refuses only
`7d75f214584a`, whose reference fit returns a boundary 0 before the gate is reached. The other
six sit at 0–24% collapse, under the quarter, so the gate passes them — including
`e7517df121d5`. **It is a guard against a reference population that is genuinely too thin, not a
discriminator among these seven**, and `e7517df121d5`'s 0.0517 is pinned with a factor-of-3.7
uncertainty that nothing downstream allows for (SEARCH_ALGORITHM §15). The earlier claim that it
"rejects exactly `e7517df121d5`" was wrong, and the next paragraph is why.

**A second cause: the check was itself random.** The gate drew its subsamples from the ambient
RNG, so it was a *coin flip re-rolled at every call* — the same reference population, unchanged,
could pin 0.0517 on one iteration and fall back to the fitted value on the next, which is
precisely the thrashing §7 exists to stop. It also meant the conclusion above was read off a
single lucky draw: repeated eight times on the full store, that build refused `e7517df121d5`
once. The draws are now seeded from a hash of the reference configurations and the caller's RNG
stream is restored on exit, so the verdict is a function of the evidence and changes only when
the evidence does. Separately, the subsample size was `max(reference_min_configs, 2n/3)`, which
with 24 configurations and a floor of 20 could drop only four — the probe was weaker than the
one the threshold was calibrated against. It is now `max(12, ceil(2n/3))`, independent of the
acceptance floor.

**Generalise.** Four compounding things, each reasonable alone. An estimate was promoted from
*a number the fit reports* to *a number the search is told to believe*, and the promotion did not
come with a check that the number could bear the weight. The sufficiency floor counted the wrong
unit — levels rather than observations within level — which is the standard way a
variance-component estimate gets trusted too early. The validation that caught it was a
resample, not a formula: when the question is "does this estimate depend on which data I
happened to get", perturbing the data answers it directly and no summary statistic does. And
then the resample was installed as a *gate* while staying random, which turns a validation into
a source of variance: **a stochastic accept/reject re-rolled on unchanged evidence is worse than
no gate**, because the quantity it guards now flickers. Seed a validation from the evidence it
is validating.

**Lives in.** `R/optimizer.R::.reference_sd_config`, `settings.R::reference_stability_draws`,
SEARCH_ALGORITHM §11.

### 45. "A new run cannot see the old rows" was false, and it had shaped the design

**Symptom.** `SEARCH_ALGORITHM` §13 stated that changing a run-defining setting "starts a fresh
run against an empty slice, leaving the old rows in the archive, unreachable from the new run",
and the whole argument about `run_id` being too broad rested on that cost. It is not true.

**What is actually true.** Selection composes exactly three filters —
`filter_evals_to_universe`, `filter_evals_to_scheme`, `filter_evals_to_build` — and **none of
them consults `run_id`**. The store is one archive shared by every run, so a new run over the
same domain, scheme and build reads every row the old run wrote. No migration and no new store
is involved in changing the run signature.

A changed `run_id` costs two real but different things: the trial universe is **re-resolved**
against the live catalogue by `pin_trial_universe`, so if a trial has been renamed or become
infeasible the new universe is smaller and `filter_evals_to_universe` *does* then gate rows out
(this is the replicate-livelock failure of #32 in another guise); and the run record splits, so
any analysis grouping by `run_id` sees one target's history as two runs.

**Why the false version was plausible.** `run_id` is stamped on every row, there is a `runs`
table keyed on it, `explain_run_mismatch()` exists to explain when it does not match, and the
tooling groups by it everywhere. Everything about the surface says it is the unit of evidence.
It is the unit of *provenance*, and nothing reads it when deciding what to learn from.

**Generalise.** A column that appears on every row, has its own table and its own diagnostic
function will be assumed to be load-bearing. Before designing around the cost of changing an
identifier, grep for what reads it — the answer here was "the report and the dev scripts", which
changes the decision completely. The fix was two settings' worth of code because the hard part
had never been hard.

**Lives in.** `R/store.R::.run_signature_keys`/`run_id_for`, `settings.R::run_identity`,
SEARCH_ALGORITHM §13.
