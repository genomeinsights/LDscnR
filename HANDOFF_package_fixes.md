# Handoff — five fixes in the LDscnR package

**Repo: `~/gitlab/LDscnR` (the R package).** This is *not* `~/gitlab/LDscnR-paper`.
Do not edit the paper repo from this session; nothing here touches the manuscript.

Five independent fixes, ordered by severity. Apply 1–4 first, verify, commit. Do 5
last and separately, because it is the only one that changes a hot numerical path.

## Ground rules

- One commit per fix. Do not batch them — fix 5 must be revertable on its own.
- After every fix: `devtools::document()` then `devtools::test()`. Roxygen is the
  source of truth for `man/` and `NAMESPACE`; never hand-edit those.
- Follow the repo's existing style: `data.table` idioms, `get("col")` inside `[`,
  `stats::`/`utils::` prefixes on base-package calls.
- If a fix turns out to be wrong or already applied, say so and stop. Do not
  invent a substitute.

---

## Fix 1 — `ld_cscore()` crashes on NA p-values

**File:** `R/ld_cscore.R`, in the `qstar` loop.

**Current:**

```r
      qv <- stats::p.adjust(p[cand], "BH")
      for (al in alpha) {
        hit <- cand[qv < al]
        if (length(hit)) cnt[hit] <- cnt[hit] + 1L
      }
```

**Problem.** `p.adjust(..., "BH")` propagates `NA`, so `qv < al` yields `NA`,
`cand[NA]` puts `NA` into `hit`, and `cnt[hit] <- ...` raises
`NAs are not allowed in subscripted assignments`. Verified:

```
hit <- cand[qv < 0.05]          #> 1, NA
cnt[hit] <- cnt[hit] + 1L       #> Error: NAs are not allowed in subscripted assignments
```

NA p-values are routine (monomorphic markers, EMMAX failures, LFMM missingness),
and this surfaces inside `ld_null_from_p()`'s surrogate loop — so it fails at
surrogate *b* of 200 after hours of compute, and under `mclapply()` it returns a
try-error rather than stopping.

**Change:** replace `cand[qv < al]` with `cand[which(qv < al)]`.

`which()` drops `NA`. This is the idiom the package already uses in
`.region_table()` in `R/ld_region_scan.R`, where the comment reads
`# which() also drops NAs`. Keep it consistent with that.

**Do not** add `na.rm`, filter NAs out of `p` earlier, or change `ncell`. A
marker with an NA p-value should simply never be a hit; its C-score stays 0 and
the denominator is unchanged.

**Verify:** add to `tests/testthat/test-ld_scan.R`:

```r
test_that("ld_cscore tolerates NA p-values", {
  set.seed(1)
  n <- 200L
  ld_ws <- matrix(runif(n * 3), n, 3,
                  dimnames = list(paste0("m", seq_len(n)), c("0.5", "0.7", "0.9")))
  p <- runif(n); p[c(3L, 17L, 200L)] <- NA_real_
  C <- expect_silent(ld_cscore(p, ld_ws))
  expect_length(C, n)
  expect_true(all(C[c(3L, 17L, 200L)] == 0))
  expect_true(all(C >= 0 & C <= 1))
})
```

---

## Fix 2 — `rho_ld` default is inconsistent across four places

**ASK THE USER BEFORE CHANGING ANYTHING HERE. This is a scientific default, not
a bug, and guessing produces different regions from the same data.**

The four disagreeing sources:

| Where | Value |
|---|---|
| `R/ld_scan.R` formal default | `rho_ld = 0.9` |
| `R/ld_scan.R` `@examples` block | `rho_ld = 0.6` |
| `R/ld_regions.R` docs | "`rho_ld = 0.9` is the canonical value" |
| every analysis in `LDscnR-paper` | `RHO_LD <- 0.60` |

Also `rho_ld = 0.9` in `R/ld_outlier_regions.R:62`, `R/ld_region_stability.R:42`,
`R/ld_regions.R:36`.

Put this to the user as two options:

- **(a) Align on 0.6** — matches every published analysis. The package is at
  `0.0.0.9000` and unreleased, so the cost is low. Change the default in all four
  R files *and* the `ld_regions.R` doc sentence, and add a `NEWS.md` entry.
- **(b) Keep 0.9** — then fix the `@examples` block in `ld_scan.R` to use 0.9, so
  the example and the default agree.

Whichever they pick, the requirement is that **all five locations state the same
value**. Do not change some and leave others.

---

## Fix 3 — no name-based alignment check

**Files:** `R/ld_cscore.R` and `R/ld_null_from_p.R`.

Both check `length(p) == nrow(ld_ws)` and never check names. A p-vector of the
right length in the wrong order gives a silently wrong answer — no error, no
warning, wrong regions.

**In `ld_cscore()`,** after the existing `stopifnot(length(p) == nrow(ld_ws))`:

```r
  if (!is.null(names(p)) && !is.null(rownames(ld_ws)) &&
      !identical(names(p), rownames(ld_ws)))
    stop("`p` is named but its names do not match `rownames(ld_ws)` -- ",
         "reorder it to the rows of `ld_ws` before calling.")
```

**In `ld_null_from_p()`,** apply the same check to `p_obs` (next to the existing
length check), and inside the `one(b)` surrogate closure to `pv`.

**Deliberately conservative:** stop, do not silently reorder. Silent reordering
would hide an upstream bug in the caller's own pipeline. Only check when names
are present — unnamed vectors stay legal, as they are now.

---

## Fix 4 — `ld_gate()` passes when the observed data has no regions

**File:** `R/ld_region_scan.R`, in `ld_gate()`.

**Current:**

```r
  out[, "ratio" := ifelse(get("obs_regions") > 0, get("med_regions") / get("obs_regions"), NA_real_)]
  out[, "pass" := is.na(get("ratio")) | get("ratio") < warn_at]
```

**Problem.** With `obs_regions == 0`, `ratio` is `NA` and `pass` becomes `TRUE`.
"Observed found nothing, surrogates found plenty" is the most emphatic gate
failure there is, and it reports as a pass.

**Change:** a basis passes only if the observed data produced regions *and* the
median surrogate stays below `warn_at` of that count:

```r
  out[, "pass" := get("obs_regions") > 0 & !is.na(get("ratio")) & get("ratio") < warn_at]
```

Leave `ratio` as `NA` — it is genuinely undefined at zero observed regions, and
the `obs_regions` column already tells the reader why `pass` is `FALSE`.

Then extend the warning so the zero-region case is not reported as a ratio
problem. Split the message: bases failing on `obs_regions == 0` get "no observed
regions at this (tau, l_min)"; bases failing on the ratio keep the existing text.

**Verify:** add a test constructing an `ld_null` whose `C_obs` is all zeros and
whose surrogates carry regions; assert `pass` is `FALSE`.

---

## Fix 5 — the BH inner loop (do this last, in its own commit)

**File:** `R/ld_cscore.R`. This is the pipeline's hot path: `ld_cscore()` runs
`B + 1` times and does `length(rho) * length(qstar)` BH computations each — at
10 x 20 x 201 that is roughly 40,000 `p.adjust()` calls over up to 790k values.

Only the *rejection set* is needed, never the adjusted-p vector.

**Replacement for the body of the `qstar` loop** (this supersedes Fix 1 in the
same lines — apply Fix 1 first anyway, so it is committed independently and the
regression test exists before you touch the arithmetic):

```r
      pc <- p[cand]
      ok <- !is.na(pc)
      ps <- sort(pc[ok])
      m  <- length(ps)
      if (!m) next
      for (al in alpha) {
        k <- max(c(0L, which(ps < al * seq_len(m) / m)))
        if (k) {
          hit <- cand[ok][pc[ok] <= ps[k]]
          cnt[hit] <- cnt[hit] + 1L
        }
      }
```

**Two details that make this exactly equivalent — get both right or the fix is
silently wrong:**

1. **Strict `<`, not `<=`, in the step-up rule.** The existing code tests
   `qv < al`. BH-adjusted `q_(i) < alpha` holds iff there exists `j >= i` with
   `p_(j) < (j/m) * alpha`, so the largest such `j` must be found with strict
   `<`. Using `<=` would reject a boundary case the current code does not.

2. **`m` counts non-NA candidates only.** `p.adjust()`'s `n = length(p)` default
   is lazily evaluated *after* it drops NAs internally, so its BH denominator is
   the number of non-NA values — which is why `m <- length(ps)` is taken from the
   NA-filtered, sorted vector. Do not use `length(cand)`.

The `pc[ok] <= ps[k]` selection (rather than taking the first `k` in sorted
order) is deliberate: it rejects all ties at the boundary together, which is what
`p.adjust()` does, since equal p-values receive equal q.

**Verify with a randomised equivalence test — do not merge without it:**

```r
test_that("the BH rejection rule matches p.adjust exactly", {
  set.seed(42)
  for (i in 1:200) {
    m  <- sample(2:300, 1L)
    p  <- c(runif(m), rep(NA_real_, sample(0:5, 1L)))
    p  <- sample(p)
    al <- sample(c(0.01, 0.05, 0.1), 1L)
    old <- which(stats::p.adjust(p, "BH") < al)
    ok  <- !is.na(p); ps <- sort(p[ok]); mm <- length(ps)
    k   <- max(c(0L, which(ps < al * seq_len(mm) / mm)))
    new <- if (k) which(ok)[p[ok] <= ps[k]] else integer(0)
    expect_identical(sort(new), sort(old))
  }
})
```

**Benchmark before and after** on a realistic shape (say 200k SNPs x 10 rho x 20
qstar) and report the speedup in the commit message. If it is under ~1.5x, say so
rather than claiming a win.

**OPTIONAL, only if the benchmark justifies it:** for a fixed `rho`, the candidate
sets are *nested* as `q*` increases, so the sort can be hoisted out of the `qstar`
loop — order `p` once per `rho` and subset that ordering. Bigger win, more room
for error. Treat it as a separate change with the same equivalence test, and skip
it if the simple version is already fast enough.

---

## Also worth doing while you are in here

`R/ld_region_c2.R`: the anchor loop is row-at-a-time with a full filter over
`sig_all` per anchor region — O(regions x significant-cells). `foverlaps()` is
already imported and used in `.region_pq()` two functions up. Vectorising it
would be consistent with the rest of the file. Not a bug; do it only if fixes
1–5 are done and verified.

`R/ld_region_c2.R` docs: `anchor_tau = 0.05` is not on the default
`tau_grid = seq(0.02, 0.5, by = 0.02)`. That is deliberate and the overlap-based
`c2` handles it correctly, but it surprises readers — add a line to `@details`.

## Final verification

```r
devtools::document()
devtools::test()
devtools::check()          # or R CMD check
```

`devtools::check()` must not gain any NOTE, WARNING or ERROR that was not there
before you started. Record the before/after check summary in the last commit
message.

## Do not

- Do not touch `~/gitlab/LDscnR-paper` from this session.
- Do not hand-edit `man/` or `NAMESPACE` — regenerate with roxygen.
- Do not change `tau`, `l_min`, `alpha`, `dcap` or `qstar` defaults. Only
  `rho_ld` is in scope, and only after the user answers Fix 2.
- Do not "improve" `.region_pq()`, `ld_gate()`'s statistics, or the C-score
  definition. The only statistical change sanctioned here is Fix 5, which must be
  provably identical to the current behaviour.
