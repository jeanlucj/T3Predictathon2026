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
| §1 The T3 data layer | #1–#10 |
| §2 Genomic methods | #11–#14 |
| §3 Search space and design | #15–#17 |
| §4 Testing and validation | #18–#20 |

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
**Now.** Recorded, not yet fixed. The remedy is replication and pairing across trials, not precision;
see the not-implemented note in `BACKGROUND.md` §4. Note also that `optimizer_step()` samples a fresh
trial per evaluation, so the trial effect is currently not even *estimable*.

### 20. Your analysis scripts deserve the same suspicion as the pipeline

**Trap.** Trusting a throwaway measurement script because the code under test is the "real" code.
**Symptom.** A comparison of two optimizer arms reported differences of *exactly* zero across nine
cells. The cause was a duplicated `n_random_init` key in a `list()` passed to `modifyList`, which
takes the first match — so both arms ran the same configuration. The tell was not an error but an
implausibly clean result.
**Now.** The comparison script asserts its own arms differ before running. The general rule is the
one behind every entry here: **trust agreement between two independent derivations; distrust a
single number** — including one your own analysis produced.
