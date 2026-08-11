# Lessons — traps we fell into, and what stops us falling in again

Every entry here is something that **looked correct, ran without error, and produced a plausible
wrong answer**. They are recorded away from the code because a fresh reader does not need our
history to understand what the code does — but *we* need it, because most of these are re-enterable
by a well-meaning simplification.

The code cites these by number (`see LESSONS.md #11`). Numbers are stable: append, never renumber.

**This is not the debugging guide.** When something looks wrong *right now*, go to `EVALUATION.md`
§6 (subtle-bug catalogue) and §9 (tooling). This file is why those checks exist.

| | |
|---|---|
| §1 The T3 data layer | #1–#10, #27 |
| §2 Genomic methods | #11–#14 |
| §3 Search space and design | #15–#17 |
| §4 Testing and validation | #18–#20 |
| §5 Numerical robustness | #21–#22 |
| §6 Running it: resources, concurrency, configuration | #23–#26, #28 |

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

### 10. A timing-out download must not be re-stormed

**Trap.** Retrying a failed download with full effort on every trial that needs it.
**Symptom.** One sick archive (e.g. project 8217) consumes an entire run, re-downloading on every
covering trial × scheme.
**Now.** Failures are tracked **per session** in RAM by project id: each attempt tries less hard
(fewer retries, shorter backoff), and after `settings$vcf_max_download_attempts` (default 3) the
project is skipped for the rest of the run. This is deliberately *not* an `unparseable_` verdict —
it is session-only and resets next run, so a project that was merely down is retried fresh. A
success clears the counter.

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
*parameters*). When adding to the search space, follow `EVALUATION.md` §10 in full. And when a
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
**Now.** Partly fixed. `trial_replication` makes a trial carry several configurations, so the trial
effect is estimable and `aggregate_scores()` removes it with random-effects BLUPs. The other half —
replicating *configurations* — is not built yet: `config_replication` is the next revision, and
until it lands every configuration gets exactly one evaluation, so the incumbent rests on a single
observation. See `BACKGROUND.md` §4.

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
**See also.** `EM_COMBINE_COMPARISON.md` — the sibling `Brapi_pipeline_for_selection` also
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
so the budget is set from measurement (`report_memory.R`) rather than from arithmetic.
**Lives in.** `R/pipeline.R::.dosage_thin_plan` / `.qc_markers`, `R/memory.R`, `report_memory.R`,
`monitor_memory.sh`.

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
**Now.** `should_backup_now()` does not consult `is_leader`, and throttles on the **backup file's
own mtime** rather than a per-process timestamp -- so N workers share one interval instead of each
honouring it separately, and the interval survives a restart. The final on-exit backup is likewise
every-worker: a leader-only one does nothing if worker 1 is the process that gets OOM-killed. The
report prints the backup's age and flags it stale past `2 x db_backup_minutes`, because
`write_report` is *not* leader-gated and so keeps updating while a long evaluation runs -- it is
the one artefact that can report the stall while the stall is happening.
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
**So watch it rather than cap it.** `report_timing.R` and `report_memory.R` exist for this.
Size a machine from the *maximum* per worker, not the median -- SLURM kills a job that exceeds
its allocation rather than throttling it -- and expect a job's final hours to lose whatever is
in flight, because a 33-hour evaluation cannot be finished by any amount of margin.
**Lives in.** `report_timing.R`, `report_memory.R`, `settings$dosage_total_budget_bytes`.
