# test_subtasks.R  --  Tier 1: deterministic core of the pipeline.
#
# Each test states its ORACLE (the expected answer, knowable without reading the
# implementation) in a comment, so reviewing the test is faster than auditing the
# code. See README_DEVELOP.md (test suite). All offline. Run: Rscript tests/test_subtasks.R

suppressMessages(library(tidyverse))
here::i_am("tests/test_subtasks.R")
source(here::here("settings.R"))
invisible(lapply(list.files(here::here("R"), "[.]R$", full.names = TRUE), source))

ok <- 0L; fail <- 0L
check <- function(cond, msg) {
  if (isTRUE(cond)) ok <<- ok + 1L
  else { fail <<- fail + 1L; cat("  FAIL:", msg, "\n") }
}
approx <- function(a, b, tol = 1e-6) all(abs(a - b) <= tol)

# Minimal config builders (only the fields each subtask reads). Overrides in `...`
# are given short names and auto-prefixed with the subtask, so e.g.
# kcfg("vanRaden_single", ridge = 0.01) sets kernel.ridge.
.prefixed <- function(prefix, ...) {
  ov <- list(...); if (length(ov)) names(ov) <- paste0(prefix, names(ov)); ov
}
kcfg <- function(method, ...) modifyList(list(
  kernel.method = method, kernel.ridge = 0, kernel.maf = 0.05,
  kernel.max_missing = 0.5, kernel.impute = "mean_round", kernel.rkhs_theta = 1),
  .prefixed("kernel.", ...))
pcfg <- function(method, ...) modifyList(list(
  pheno_prep.method = method, pheno_prep.z_thr = Inf, pheno_prep.standardize = "none",
  pheno_prep.ge_weighting = "none", pheno_prep.ge_bandwidth = 1),
  .prefixed("pheno_prep.", ...))
prcfg <- function(method, ...) modifyList(list(
  predict_post.method = method, predict_post.blend_obs_w = 0, predict_post.min_overlap = 2),
  .prefixed("predict_post.", ...))
no_trial <- list(lat = NA_real_, long = NA_real_, year = NA_integer_)

# ===========================================================================
cat(".brapi_try (retry through a flaky server)\n")
# Oracle: retries a transient error and returns the eventual success; a persistent error
# raises after exactly `tries` attempts; tries=1 disables retry. base_delay=0 keeps it fast.
n <- 0L; th <- function() { n <<- n + 1L; if (n < 3L) stop("boom"); "ok" }
check(identical(suppressMessages(.brapi_try(th, tries = 4, base_delay = 0)), "ok") && n == 3L,
      "retries a transient error, then succeeds (3rd attempt)")
n2 <- 0L; th2 <- function() { n2 <<- n2 + 1L; stop("always") }
check(inherits(suppressMessages(try(.brapi_try(th2, tries = 3, base_delay = 0), silent = TRUE)),
               "try-error") && n2 == 3L, "persistent error raises after exactly `tries` attempts")
n3 <- 0L; th3 <- function() { n3 <<- n3 + 1L; stop("x") }
suppressMessages(try(.brapi_try(th3, tries = 1, base_delay = 0), silent = TRUE))
check(n3 == 1L, "tries=1 disables retry (single attempt)")

# ===========================================================================
cat("score_predictions\n")
# Oracle: correlation of a vector with a positive affine transform of itself is 1.
obs <- setNames(c(1, 2, 3, 4, 5, 6), letters[1:6])
r <- score_predictions(setNames(2 * obs + 10, names(obs)), obs)
check(r$status == "ok" && r$n_test == 6 && approx(r$score, 1), "self-correlation = 1")
# Oracle: a negative transform gives -1.
r <- score_predictions(setNames(-obs, names(obs)), obs)
check(approx(r$score, -1), "negated = -1")
# Oracle: only shared accession names are joined (extras on each side ignored).
p <- c(a = 1, b = 2, c = 3, d = 4, e = 5, f = 6, z = 99)   # z not in obs
o <- c(obs, y = 99)                                         # y not in pred
r <- score_predictions(p, o)
check(r$n_test == 6, "join keeps only shared names")
# Oracle: fewer than 5 overlapping accessions is a defined failure.
r <- score_predictions(c(a = 1, b = 2, c = 3, d = 4), obs)
check(is.na(r$score) && r$status == "too_few_overlap", "n<5 -> too_few_overlap")
# Oracle: constant predictions have no variance -> defined failure.
r <- score_predictions(setNames(rep(5, 6), names(obs)), obs)
check(is.na(r$score) && r$status == "constant", "constant pred -> constant")

# ===========================================================================
cat(".vcf_to_dosage (orientation + encoding; guards the transposition bug)\n")
vcf_lines <- c(
  "##fileformat=VCFv4.2",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
        "S1", "S2", sep = "\t"),
  paste("1", "100", "m1", "A", "G", ".", ".", ".", "GT", "0/0", "1/1", sep = "\t"),
  paste("1", "200", "m2", "A", "G", ".", ".", ".", "GT", "0/1", "0/0", sep = "\t"),
  paste("1", "300", "m3", "A", "G", ".", ".", ".", "GT", "1/1", "0/1", sep = "\t"))
vpath <- tempfile(fileext = ".vcf"); writeLines(vcf_lines, vpath)
# Oracle (hand-decoded): S1 = (m1=0/0=0, m2=0/1=1, m3=1/1=2); S2 = (2, 0, 1).
d <- suppressWarnings(.vcf_to_dosage(vpath, c("S1", "S2"), thin = 1))
check(identical(dim(d), c(2L, 3L)), "dosage is accessions x markers (2x3)")
check(identical(rownames(d), c("S1", "S2")), "rownames are samples")
check(identical(colnames(d), c("m1", "m2", "m3")), "colnames are marker IDs")
check(approx(d["S1", ], c(0, 1, 2)) && approx(d["S2", ], c(2, 0, 1)), "encoding + orientation")
# Oracle: subsetting to one sample yields one row.
d1 <- suppressWarnings(.vcf_to_dosage(vpath, "S2", thin = 1))
check(nrow(d1) == 1 && approx(d1["S2", ], c(2, 0, 1)), "single-sample subset")
# Oracle: thin=2 keeps markers 1 and 3 (every 2nd, genome-wide).
dt <- suppressWarnings(.vcf_to_dosage(vpath, c("S1", "S2"), thin = 2))
check(identical(colnames(dt), c("m1", "m3")), "marker thinning keeps every 2nd")

# Oracle: a malformed variant line (wrong field count) is SKIPPED, not fatal -- the good
# markers still come through. (This is the real 8114 failure: one 3-field variant record.)
bad_lines <- c(vcf_lines[1:5],
               paste("1", "250", "BAD", sep = "\t"),            # 3 fields, not 11
               vcf_lines[6])
bpath <- tempfile(fileext = ".vcf"); writeLines(bad_lines, bpath)
db <- .vcf_to_dosage(bpath, c("S1", "S2"), thin = 1)
check(identical(colnames(db), c("m1", "m2", "m3")), "malformed variant line skipped, good markers kept")

# Oracle: a transposed / non-VCF header (markers as columns) is REJECTED, not misread.
# (The real 14511: the "#CHROM" row holds chromosome labels; field 2 is not POS.)
tr_lines <- c("##fileformat=VCFv4.2",
              paste("#CHROM", "1A", "1A", "1B", "2A", sep = "\t"),
              paste("POS", "100", "200", "300", "400", sep = "\t"))
tpath <- tempfile(fileext = ".vcf"); writeLines(tr_lines, tpath)
check(inherits(try(.vcf_to_dosage(tpath, NULL, 1L), silent = TRUE), "try-error"),
      "transposed / non-VCF header rejected")

# Oracle: .vcf_stat validates + counts in one pass.
st <- .vcf_stat(vpath)
check(st$n_samples == 2 && st$n_markers == 3 && identical(st$samples, c("S1", "S2")),
      ".vcf_stat: 2 samples, 3 markers")
check(inherits(try(.vcf_stat(tpath), silent = TRUE), "try-error"), ".vcf_stat rejects non-VCF")

# Oracle: .eff_thin -- fits under budget -> requested thin; over budget -> raised to fit.
# Integer inputs on purpose: n_samples*n_markers must not overflow 32-bit int (683L*7.47ML).
check(.eff_thin(100L, 1000L, 1L, 2e9) == 1, ".eff_thin: small project keeps requested thin")
check(.eff_thin(100L, 1000L, 4L, 2e9) == 4, ".eff_thin: never below the requested thin")
# 683 x 7.47M x 4 = 20.4 GB; /2 GB budget -> ceil 10.2 -> 11.
check(.eff_thin(683L, 7467224L, 1L, 2e9) == 11, ".eff_thin: oversized project auto-thinned (no int overflow)")

# ===========================================================================
cat(".qc_markers\n")
set.seed(1)
good <- matrix(rep(c(0, 1, 2, 0, 1, 2), 60), nrow = 6)        # 60 polymorphic markers
colnames(good) <- paste0("g", 1:60)
Xqc <- cbind(good,
             mono = rep(0, 6),                  # MAF 0 -> dropped
             hm   = c(0, 1, NA, NA, NA, NA),     # 4/6 missing > 0.5 -> dropped
             imp  = c(0, 2, 2, 2, 2, NA))        # one NA -> imputed
rownames(Xqc) <- paste0("id", 1:6)
q <- .qc_markers(Xqc, kcfg("vanRaden_single"))
check(!("mono" %in% colnames(q)), "monomorphic marker dropped")
check(!("hm" %in% colnames(q)), "high-missing marker dropped")
check("imp" %in% colnames(q) && !anyNA(q), "kept marker imputed (no NA remains)")
# Oracle: mean of non-NA imp values is (0+2+2+2+2)/5 = 1.6; mean_round -> 2.
check(q["id6", "imp"] == 2, "mean_round imputation = round(1.6) = 2")
check(.qc_markers(Xqc, kcfg("vanRaden_single", impute = "mean"))["id6", "imp"] == 1.6,
      "mean imputation = 1.6")
# Oracle: fewer than 50 surviving markers is an error.
check(inherits(try(.qc_markers(good[, 1:10], kcfg("vanRaden_single")), silent = TRUE),
               "try-error"), "errors below 50 markers")

# ===========================================================================
cat(".merge_markers\n")
mk <- function(rows, cols) { m <- matrix(1, length(rows), length(cols),
  dimnames = list(rows, cols)); m }
d1 <- mk(c("a", "b", "c", "d"), paste0("m", 1:60))
d2 <- mk(c("c", "d", "e", "f"), paste0("m", 11:80))   # shares m11..m60 = 50 markers
mg <- .merge_markers(list(d1, d2))
# Oracle: 50 shared markers; accessions a,b,c,d,e,f with duplicates c,d dropped.
check(ncol(mg) == 50, "merge keeps the 50-marker intersection")
check(setequal(rownames(mg), c("a", "b", "c", "d", "e", "f")) && !anyDuplicated(rownames(mg)),
      "every accession once; duplicates c,d dropped")
# Oracle: when an accession sits in two projects, its row must come from the project with
# the RICHER marker build (more markers = the better-genotyped call of the same line).
rich <- matrix(2, 2, 70, dimnames = list(c("c", "d"), paste0("m", 11:80)))   # 70 markers, all 2s
poor <- matrix(0, 2, 60, dimnames = list(c("c", "d"), paste0("m", 1:60)))    # 60 markers, all 0s
check(all(.merge_markers(list(poor = poor, rich = rich))["c", ] == 2),
      "duplicated accession takes its row from the richer marker build")
# Oracle: <50 shared markers -> fall back to the project with more accessions.
small1 <- mk(c("a", "b", "c", "d"), paste0("m", 1:40))
small2 <- mk(c("a", "b", "c", "d", "e", "f"), paste0("p", 1:40))   # 0 shared markers
check(identical(.merge_markers(list(small1, small2)), small2), "no-overlap -> largest project")
check(identical(.merge_markers(list(d1)), d1), "single project returned as-is")

# ===========================================================================
cat("build_kernel\n")
set.seed(2)
Xc <- matrix(sample(0:2, 6 * 60, replace = TRUE), nrow = 6)
Xc[2, ] <- Xc[1, ]                                    # id2 is a clone of id1
rownames(Xc) <- paste0("id", 1:6); colnames(Xc) <- paste0("m", 1:60)
all6 <- rownames(Xc)
K <- build_kernel(kcfg("vanRaden_single"), list(Xc), all6)
check(isSymmetric(unname(K)), "VanRaden GRM symmetric")
# Oracle: genetically identical lines have identical relationship vectors.
check(approx(K["id1", ], K["id2", ], tol = 1e-8), "clone lines -> identical GRM rows")
# Oracle: adding ridge raises only the diagonal, by exactly ridge.
K0 <- build_kernel(kcfg("vanRaden_single", ridge = 0), list(Xc), all6)
K1 <- build_kernel(kcfg("vanRaden_single", ridge = 0.01), list(Xc), all6)
check(approx(diag(K1) - diag(K0), 0.01, 1e-8), "ridge adds exactly 0.01 to diagonal")
check(approx((K1 - K0)[upper.tri(K1)], 0, 1e-8), "ridge leaves off-diagonals unchanged")
# Oracle: em_combine with a single project is the plain VanRaden GRM.
check(approx(build_kernel(kcfg("em_combine"), list(Xc), all6), K0, 1e-8),
      "em_combine, one project -> VanRaden")
# RKHS Gaussian kernel.
Kr <- build_kernel(kcfg("rkhs_gaussian", rkhs_theta = 1, ridge = 0), list(Xc), all6)
check(isSymmetric(unname(Kr)), "RKHS kernel symmetric")
check(approx(diag(Kr), 1, 1e-8), "RKHS diagonal = 1")
check(all(Kr > 0 & Kr <= 1 + 1e-9), "RKHS entries in (0,1]")
check(approx(Kr["id1", "id2"], 1, 1e-8), "RKHS clone pair = 1")

# ===========================================================================
cat(".vanraden (population frequencies; the A.mat coding bug)\n")
set.seed(7)
npop <- 200; nmk <- 300
Xp <- matrix(sample(0:2, npop * nmk, replace = TRUE), nrow = npop,
             dimnames = list(paste0("g", 1:npop), paste0("m", 1:nmk)))
need4 <- c("g3", "g17", "g42", "g99")
# Oracle (the whole point of the change): the GRM built over `need` with population
# allele frequencies IS the [need, need] submatrix of the whole population's GRM.
# So "QC + GRM on the population, then subset" and this are the same numbers.
Kfull <- .vanraden(Xp, rownames(Xp))
Ksub  <- .vanraden(Xp, need4)
check(approx(Ksub, Kfull[need4, need4], 1e-10),
      "GRM on needed rows w/ population freqs == submatrix of the full-population GRM")

# Oracle: hand-computed VanRaden on a panel where EVERY marker is alt-major (p > 0.5).
# This is the case rrBLUP::A.mat silently destroys: it assumes {-1,0,1}, so on {0,1,2}
# its internal freq = p + 0.5 > 1, its MAF goes negative, and min.MAF drops the marker.
set.seed(9)
pmaj <- runif(nmk, 0.7, 0.95)                        # every marker alt-major
Xmaj <- sapply(pmaj, function(pk) rbinom(60, 2, pk))
dimnames(Xmaj) <- list(paste0("g", 1:60), paste0("m", 1:nmk))
phat  <- colMeans(Xmaj) / 2
Wman  <- sweep(Xmaj, 2, 2 * phat, "-")
Kman  <- tcrossprod(Wman) / (2 * sum(phat * (1 - phat)))   # VanRaden, by hand
check(approx(.vanraden(Xmaj, rownames(Xmaj)), Kman, 1e-10),
      "alt-major panel: .vanraden matches hand-computed VanRaden")
# ... and the regression itself: A.mat, fed the same {0,1,2} matrix, does NOT. It reads
# every alt-major marker as MAF < 0 and drops it, so on this panel it drops ALL of them
# and returns an entirely NaN GRM. (On real panels it drops ~half and returns a GRM that
# merely looks plausible -- which is why this went unnoticed.)
Abad <- suppressWarnings(tryCatch(rrBLUP::A.mat(Xmaj), error = function(e) NULL))
check(is.null(Abad) || anyNA(Abad) || !isTRUE(approx(Abad, Kman, 1e-6)),
      "regression: rrBLUP::A.mat on {0,1,2} dosages does not give VanRaden (the bug)")
check(is.null(Abad) || all(is.na(Abad)),
      "regression: A.mat drops every alt-major marker -> all-NaN GRM")

# Oracle: the {0,1,2} assumption now lives in .vanraden's `p`, so it must be guarded.
# Handing it rrBLUP's {-1,0,1} coding must FAIL LOUDLY, not silently return a wrong p.
check(inherits(try(.vanraden(Xp - 1, rownames(Xp)), silent = TRUE), "try-error"),
      "guard: {-1,0,1} coding is rejected, not silently mis-centred")

# Oracle: QC decisions depend on the population they are estimated from. QC on a
# 5-accession subset keeps a different marker set than QC on the population.
qc_pop <- colnames(.qc_markers(Xp, kcfg("vanRaden_single")))
qc_sub <- colnames(.qc_markers(Xp[1:5, , drop = FALSE], kcfg("vanRaden_single")))
check(!identical(qc_pop, qc_sub), "QC on 5 accessions != QC on the 200-accession population")

# ===========================================================================
cat(".bridge_accessions + em_combine stitching\n")
set.seed(11)
mkpanel <- function(rows, mcols) {
  m <- matrix(sample(0:2, length(rows) * length(mcols), replace = TRUE),
              nrow = length(rows), dimnames = list(rows, mcols))
  m
}
# Two DIFFERENT protocols (disjoint marker sets). need = {n1, n2}: n1 only in panel A,
# n2 only in panel B. b1..b3 are genotyped in BOTH -> the bridge.
bridge3 <- paste0("b", 1:3)
pA <- mkpanel(c("n1", bridge3, paste0("xa", 1:8)), paste0("mA", 1:120))   # protocol A markers
pB <- mkpanel(c("n2", bridge3, paste0("xb", 1:8)), paste0("mB", 1:120))   # protocol B markers
need2 <- c("n1", "n2")

# Oracle: bridge = accessions in >= 2 panels; disjoint panels -> none.
check(setequal(.bridge_accessions(list(pA, pB)), bridge3),
      ".bridge_accessions finds accessions in >= 2 panels")
check(length(.bridge_accessions(list(pA, mkpanel(paste0("z", 1:6), paste0("mC", 1:120))))) == 0,
      "disjoint panels -> no bridge accessions")

# Oracle (the crux): with need alone, panel A has a single needed row (n1) -> .vanraden NULL,
# uncombinable. Keeping the bridge accessions makes the panel a valid GRM.
check(is.null(.vanraden(.qc_markers(pA, kcfg("em_combine")), need2)),
      "need-only: panel with one needed accession is uncombinable (NULL)")
check(!is.null(.vanraden(.qc_markers(pA, kcfg("em_combine")), union(need2, bridge3))),
      "with bridges: the same panel yields a valid GRM")

# Oracle: em_combine end to end -- feasible BECAUSE of the bridges, and the bridges are
# eliminated so the result is exactly need x need. (Runs the real covariance_combiner.)
Kem <- build_kernel(kcfg("em_combine", ridge = 0), list(A = pA, B = pB), need2)
check(setequal(rownames(Kem), need2), "em_combine result is exactly need x need (bridges dropped)")
check(isSymmetric(unname(Kem)), "em_combine stitched GRM is symmetric")
check(all(is.finite(Kem)), "em_combine stitched GRM is finite (cross-panel block estimated)")

# ===========================================================================
cat(".group_by_panel / .prune_redundant\n")
gset <- optimizer_settings()
pmk <- function(rows, cols) matrix(1, length(rows), length(cols),
                                   dimnames = list(rows, cols))
# Same panel (v1 vs v2 reference build): v2 contains ALL of v1's markers.
v1 <- pmk(paste0("a", 1:10), paste0("m", 1:100))
v2 <- pmk(paste0("a", 1:10), paste0("m", 1:140))       # 100/100 = 100% containment
# A different panel: shares only 40 of the smaller 100 markers.
other <- pmk(paste0("b", 1:10), c(paste0("m", 1:40), paste0("x", 1:60)))
g <- .group_by_panel(list(v1 = v1, v2 = v2, other = other), gset)
check(any(vapply(g, function(x) setequal(x, c("v1", "v2")), logical(1))),
      "100% marker containment -> same protocol group")
check(any(vapply(g, function(x) identical(x, "other"), logical(1))),
      "40% marker containment -> its own group")

# Oracle: identical accession sets -> keep the richer marker build (the v1/v2 case).
pr <- .prune_redundant(list(v1 = v1, v2 = v2), gset)
check(identical(names(pr), "v2"), "identical accessions -> project with more markers survives")
# Oracle: >=90% shared accessions but not identical -> keep the one with MORE accessions,
# even though it has fewer markers.
big   <- pmk(paste0("a", 1:100), paste0("m", 1:100))
small <- pmk(paste0("a", 1:95),  paste0("m", 1:140))   # 95/95 = 100% of the smaller set
pr2 <- .prune_redundant(list(big = big, small = small), gset)
check(identical(names(pr2), "big"), "92% shared accessions -> project with more accessions survives")
# Oracle: distinct accession sets are NOT redundant.
indep <- pmk(paste0("z", 1:50), paste0("m", 1:100))
check(length(.prune_redundant(list(big = big, indep = indep), gset)) == 2,
      "disjoint accessions -> both projects kept")

# ===========================================================================
cat("build_targets + .blue_per_trial\n")
mkobs <- function(study, germ, value) tibble::tibble(
  study_id = study, germplasm_name = germ, value = value, rep = "r1", block = "b1")
# raw_mean. Oracle: mean of per-trial means. g1: trial A mean(10,12)=11, trial B=20 -> 15.5; g2=8.
obs_rm <- dplyr::bind_rows(
  mkobs("A", "g1", 10), mkobs("A", "g1", 12), mkobs("A", "g2", 8), mkobs("B", "g1", 20))
t_rm <- build_targets(pcfg("raw_mean"), obs_rm, no_trial, NULL, list())
check(approx(t_rm["g1"], 15.5) && approx(t_rm["g2"], 8), "raw_mean = mean of per-trial means")
# trial_center. Oracle: removing trial means leaves only a global shift, so adding
# a constant to one trial shifts ALL targets equally (variance of the change = 0).
t1 <- build_targets(pcfg("trial_center"), obs_rm, no_trial, NULL, list())
obs_shift <- obs_rm |> dplyr::mutate(value = value + ifelse(study_id == "B", 100, 0))
t2 <- build_targets(pcfg("trial_center"), obs_shift, no_trial, NULL, list())
common <- intersect(names(t1), names(t2))
check(approx(stats::var(t2[common] - t1[common]), 0, 1e-8), "trial_center removes trial effect")
# blue_lm. Oracle: single-level design -> germ coefficients are the group means.
obs_blue <- dplyr::bind_rows(
  mkobs("A", "g1", 10), mkobs("A", "g1", 10), mkobs("A", "g2", 20), mkobs("A", "g2", 20))
t_bl <- build_targets(pcfg("blue_lm"), obs_blue, no_trial, NULL, list())
check(approx(t_bl["g1"], 10) && approx(t_bl["g2"], 20), "blue_lm recovers group means")
# blue_lm fallback. Oracle: a one-germplasm trial -> mean. g1 = mean(10,14) = 12.
t_fb <- build_targets(pcfg("blue_lm"),
  dplyr::bind_rows(mkobs("A", "g1", 10), mkobs("A", "g1", 14)), no_trial, NULL, list())
check(approx(t_fb["g1"], 12), "blue_lm single-line fallback = mean")
# outlier removal. Oracle: with enough normal points, z_thr=3 flags the 500; g1 -> 10.
obs_out <- dplyr::bind_rows(
  mkobs("A", "g1", rep(10, 20)), mkobs("A", "g1", 500), mkobs("A", "g2", rep(20, 5)))
t_o3  <- build_targets(pcfg("blue_lm", z_thr = 3),   obs_out, no_trial, NULL, list())
t_oInf <- build_targets(pcfg("blue_lm", z_thr = Inf), obs_out, no_trial, NULL, list())
check(approx(t_o3["g1"], 10, 1e-6), "z_thr=3 removes the outlier (g1 -> 10)")
check(t_oInf["g1"] > 30, "z_thr=Inf keeps the outlier (g1 inflated)")
# two_stage_blup. Oracle: BLUP shrinks the spread below the raw germ effects.
if (requireNamespace("lme4", quietly = TRUE)) {
  obs_ts <- dplyr::bind_rows(
    mkobs("A", "g1", 4), mkobs("A", "g2", 11), mkobs("A", "g3", 14),
    mkobs("B", "g1", 6), mkobs("B", "g2", 9),  mkobs("B", "g3", 16))
  t_ts <- suppressWarnings(build_targets(pcfg("two_stage_blup"), obs_ts, no_trial, NULL, list()))
  raw  <- c(g1 = 5, g2 = 10, g3 = 15)             # per-germ means across the two trials
  check(stats::sd(t_ts - mean(t_ts)) < stats::sd(raw - mean(raw)),
        "two_stage_blup shrinks spread below raw effects")
} else cat("  (skipped two_stage_blup: lme4 not installed)\n")
# standardize. Oracle: per-trial z-scores have mean 0 and sd 1.
obs_sd <- mkobs("A", paste0("g", 1:5), 1:5)
t_sd <- build_targets(pcfg("raw_mean", standardize = "per_trial_z"), obs_sd, no_trial, NULL, list())
check(approx(mean(t_sd), 0, 1e-9) && approx(stats::sd(t_sd), 1, 1e-3), "per_trial_z -> mean 0, sd 1")

# ===========================================================================
cat("predict_test\n")
# Conditional expectation. Oracle: a test line that is a clone of training line t1
# (its K row to the training set equals t1's) gets prediction mu + u[t1].
Kp <- matrix(0, 3, 3, dimnames = list(c("t1", "t2", "x"), c("t1", "t2", "x")))
diag(Kp) <- 2; Kp["t1", "t2"] <- Kp["t2", "t1"] <- 0.5
Kp["x", "t1"] <- Kp["t1", "x"] <- 2; Kp["x", "t2"] <- Kp["t2", "x"] <- 0.5   # x clone of t1
fit <- list(kind = "gblup", u = c(t1 = 1, t2 = -0.5), mu = 10)
pc <- predict_test(prcfg("cond_expectation"), fit, Kp, c("t1", "t2"), "x",
                   targets = c(t1 = 0, t2 = 0), scheme = "CV00")
check(approx(pc["x"], 10 + 1), "cond_expectation clone -> mu + u[t1]")
# Blending. Oracle: under CV0, blend = w*BLUE + (1-w)*genomic; under CV00 it's a no-op.
Kd <- matrix(0, 4, 4, dimnames = list(c("a", "t1", "t2", "t3"), c("a", "t1", "t2", "t3")))
fitd <- list(kind = "gblup", u = c(a = 0, t1 = 0, t2 = 0, t3 = 0), mu = 10)
tg <- c(a = 100, t1 = 0, t2 = 0, t3 = 0)
p0 <- predict_test(prcfg("direct_blup", blend_obs_w = 0.5), fitd, Kd,
                   c("t1", "t2", "t3"), "a", tg, "CV0")
p00 <- predict_test(prcfg("direct_blup", blend_obs_w = 0.5), fitd, Kd,
                    c("t1", "t2", "t3"), "a", tg, "CV00")
check(approx(p0["a"], 0.5 * 100 + 0.5 * 10), "CV0 blends observed BLUE")
check(approx(p00["a"], 10), "CV00 blend is a no-op")
# Fallback. Oracle: too little genotyped overlap -> every prediction is mean(targets).
pf <- predict_test(prcfg("cond_expectation", min_overlap = 3), fit, Kp,
                   c("t1", "t2"), "x", targets = c(t1 = 4, t2 = 8), scheme = "CV0")
check(approx(pf["x"], 6), "fallback below min_overlap -> mean(targets)")

# ===========================================================================
cat("mask_cv\n")
mo <- dplyr::bind_rows(
  mkobs("A", "f1", 1), mkobs("A", "f2", 2), mkobs("A", "x1", 3), mkobs("B", "x2", 4))
focal <- c("f1", "f2")
m0  <- mask_cv(mo, focal, "CV0")
m00 <- mask_cv(mo, focal, "CV00")
check(nrow(m0) == 4, "CV0 keeps all rows")
check(nrow(m00) == 2 && !any(m00$germplasm_name %in% focal), "CV00 drops focal-accession rows")
check(setequal(setdiff(m0$germplasm_name, m00$germplasm_name), focal),
      "CV0 minus CV00 = exactly the focal accessions")

# ===========================================================================
cat(sprintf("\nTier 1 subtask tests: %d passed, %d failed\n", ok, fail))
if (fail > 0) quit(status = 1)
