# Fallback: running on SciNet's R module, no container

*Use this only if the container route stalls. The container is preferred because it pins
package versions; this route pins nothing and will drift.*

SciNet has **`r/4.5.3`**, which is why `optimizer.def` pins the container to R 4.5.3 too — the
two routes then run the same R, and switching between them cannot by itself move a result.

```bash
module spider r          # authoritative list; the public docs page lags
module load r/4.5.3
```

## The trap that will bite you here

SciNet documents installing personal packages by putting this in `.Renviron`:

```
R_LIBS_USER=/project/<account>/R_packages/%v
```

**R reads exactly one `.Renviron`** — the working directory's if there is one, otherwise
`$HOME`'s. It does not merge them and it does not walk parents.

This project *always* runs from the optimizer directory, and that directory has its own
`.Renviron` holding `T3_USERNAME` and `T3_PASSWORD`. So a `~/.Renviron` containing
`R_LIBS_USER` is **never read**, and packages land in the default library instead — on Ceres,
inside your 30 GB home.

Put both in `optimizer/.Renviron`:

```
T3_USERNAME=...
T3_PASSWORD=...
OPTIMIZER_HOME=/project/<account>/t3_optimizer
R_LIBS_USER=/project/<account>/R_packages/%v
```

Verify before installing anything:

```bash
cd <repo>/scripts/Analysis_Claude/optimizer
Rscript -e '.libPaths()'        # your project path must be FIRST
```

## Installing the packages

```bash
mkdir -p /project/<account>/R_packages/4.5
cd <repo>/scripts/Analysis_Claude/optimizer

Rscript -e 'install.packages(c(
  "tidyverse","here","DBI","RSQLite","jsonlite","rlang","Matrix",
  "lme4","rrBLUP","sommer","rpart","MASS","janitor","httr",
  "RhpcBLASctl","remotes"),
  repos = "https://packagemanager.posit.co/cran/2026-08-01")'

# Pinned to the same commits as optimizer.def -- R/data_access.R reaches a PRIVATE function
# in T3BrapiHelpers, so a branch install can break it with no warning.
Rscript -e 'remotes::install_github("TriticeaeToolbox/BrAPI.R@51d8d450d8ec4f9fc13248165b6382f4a24030b0", upgrade = "never")'
Rscript -e 'remotes::install_github("jeanlucj/T3BrapiHelpers@6c756462b5a315a992bdd7a26585d912a5452013", upgrade = "never")'
```

Expect 30–60 minutes: `lme4`, `sommer` and `Matrix` all compile.

## Then verify

```bash
cd <repo>/scripts/Analysis_Claude/optimizer
Rscript -e 'cat(Sys.which("rsync"), "\n")'   # must be non-empty; the cache backup
                                             # skips SILENTLY without it
Rscript tests/run_all.R                      # must reach 32 / 8007 / 317
Rscript peek_failures.R                      # exercises login + catalogue + store, read-only
```

`rsync` is normally present on a cluster, unlike in a bare container image — but check, because
`sync_cache_to_backup()` prints "rsync not on PATH -- skipping" and carries on with no cache
backup at all.

## What this route does not give you

The container records its R version, CRAN snapshot date and both GitHub SHAs in `%labels`, so
a score traces to an exact environment. Here, `install.packages()` at a later date silently
gives different versions. If you run this route for anything that ends up in the manuscript,
record `sessionInfo()` alongside `OPTIMIZER_BUILD`.
