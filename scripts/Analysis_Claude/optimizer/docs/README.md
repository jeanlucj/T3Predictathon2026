# `docs/` — understanding the system

Six documents, each answering a different question. Read the one that matches yours rather
than starting at the top.

| document | the question it answers |
|---|---|
| `DESIGN.md` | **How is it built?** The architecture: every function, what it returns, which subtask it belongs to, and where its outputs go. Read this before changing code. |
| `SEARCH_ALGORITHM.md` | **How does the search decide what to try next?** The full rules: every term defined, the variance components and what they mean, the contender threshold and the top-8 cap, how the budget is prioritized, what "churn" is and how to measure it — each rule anchored to the code and to what the stored runs measure it doing. §13 is where to go for **what counts as one run** -- a target you keep adding evidence to, or an experiment in one fixed policy -- and why changing that needs no store migration. Read this before concluding anything about why a run is behaving as it is. |
| `BACKGROUND.md` | **Why is it built that way?** The statistical and data-management problems — sparse marker overlap, cross-trial G×E, the CV0/CV00 distinction, why the search is only nominally better than random at realistic budgets — and how each is addressed. Read this before changing the *method*. |
| `EVALUATION.md` | **Does it work?** The evaluation runbook: exercise one module at a time (`arm_evaluation`), run and extend the tests, trace a single evaluation end to end, and use the validation tooling (the canary oracle, `diagnose_trial`). |
| `EVALUATION_CHECKLIST.md` | **Is it *correct*?** A do-list to work through before trusting a run's results. A slower question than whether it launches, and a different one from whether it runs. |
| `LESSONS.md` | **What has already gone wrong?** A numbered record of every failure that cost real time, with the diagnosis and the fix. Code comments cite it by number (`LESSONS #24`), so this is the file those pointers land in. |

`DESIGN.md` §3 and `SEARCH_ALGORITHM.md` are the summary and the full statement of the same
thing: if §3 answers your question, stop there; if you are about to change a replication knob,
read the long one first.

`DESIGN.md` versus `BACKGROUND.md` is the pair people open both of: **`DESIGN.md` is the
mechanism, `BACKGROUND.md` is the justification.** If you want to know what
`covariance_combiner` does, that is `DESIGN.md`; if you want to know why combining panels is
worth the trouble at all, that is `BACKGROUND.md`.

`LESSONS.md` is append-only and is where history goes. A pitfall gets a one-line comment in the
code plus a pointer here — not the whole story inline.

Operating a run is `../tools/README.md`; the launch do-lists are `../RUNBOOK_INTERACTIVE.md` and
`../RUNBOOK_SLURM.md`; changing the algorithm is `../dev/README.md`.
