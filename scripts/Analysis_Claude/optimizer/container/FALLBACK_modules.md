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

One script does it:

```bash
cd <repo>/scripts/Analysis_Claude/optimizer
module load r/4.5.3
Rscript setup_fallback_libs.R --dry-run     # what it would install, and from where
Rscript setup_fallback_libs.R               # 30-60 min: lme4, sommer, Matrix compile
```

It **reads the package list, the CRAN snapshot date and both git SHAs out of
`container/optimizer.def`** rather than restating them, so this route and the container cannot
drift apart. If that parse ever fails it stops rather than quietly installing today's CRAN.

It also refuses in the two situations where it would do harm: run inside the container (where
everything is already present, and anything installed lands in the ephemeral overlay), or with
`R_LIBS_USER` unset or pointing into your 30 GB home.

Afterwards it verifies with the same assertions `optimizer.def`'s `%test` makes — every package
loads, and the private `make_row_from_trial_result` still exists in `T3BrapiHelpers`.

> The script is deliberately written in **base R only**: no `library(tidyverse)`, no
> `here::i_am()`. It runs precisely when those packages are missing, so depending on them would
> make it fail in the one situation it exists for.

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
