# Cross-machine development & analysis workflow

How to develop `LDscnR` and use it in analyses across more than one computer,
without copying code or reinstalling from GitHub on every change.

## Core idea

Each machine has **one git clone** of this repo. That single clone serves both
jobs:

- **run** the package in analyses via `devtools::load_all()` (no install), and
- **edit** the package to fix bugs.

Fixes move between machines through `git push` / `git pull` — never by copying
files. Git carries the exact changes, with history, in both directions.

## One-time setup on each machine

```bash
git clone https://github.com/genomeinsights/LDscnR.git
cd LDscnR
git switch overhaul        # the branch under active development
```

In analysis scripts (e.g. in the separate LDscnR-paper project), replace
`library(LDscnR)` with:

```r
devtools::load_all("~/path/to/LDscnR")   # loads current SOURCE, no install
```

Now the code you run *is* the code you edit — bugs reproduce against real data.

## Development loop (on whichever machine reproduces the bug)

```
1. Reproduce the bug in your analysis        (load_all)
2. Edit the function in R/…
3. devtools::load_all()  →  re-test          (repeat 2–3 until fixed)
4. devtools::document()  if you changed roxygen comments
5. git commit -am "Fix <bug>"
6. git push
```

Then on the other machine:

```bash
git pull        # the fix arrives — no copying, full history
```

## Keeping the two clones from diverging

- `git pull` **before** you start editing, and again **before** you push.
- For a solo two-machine setup, working directly on `overhaul` is fine.
- For anything involved or experimental, branch it:
  `git switch -c bugfix/<name>`, push, then merge into `overhaul` when happy.

## Running vs. reproducing (load_all) — and pinning (install)

- While **bug-hunting**, `load_all()` a branch: instant, reflects every edit.
- For **reproducible paper runs**, install a *pinned tag* instead, so results
  don't shift when the branch moves on:

  ```r
  # after committing + pushing a fix, from the package repo:
  #   git tag v0.1.1 && git push origin v0.1.1
  remotes::install_github("genomeinsights/LDscnR@v0.1.1")
  ```

Current frozen release: **`v0.1.0-overhaul`** (decay core + auto MAF filtering +
`compute_ld_w`).

## Gotcha: HTTP 401 "Bad credentials" on install_github

This repo is **public**, so installing needs no token. A 401 means the machine
has a stale `GITHUB_PAT` / `GITHUB_TOKEN` being sent and rejected. Clear it:

```r
Sys.unsetenv("GITHUB_PAT"); Sys.unsetenv("GITHUB_TOKEN")
remotes::install_github("genomeinsights/LDscnR@v0.1.0-overhaul")
```

Permanent fix: remove/replace the `GITHUB_PAT` line via
`usethis::edit_r_environ()`, then restart R.
