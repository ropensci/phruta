# Changelog

## `phruta` 0.1.4

- Fixed a bug in
  [`sq.retrieve.direct()`](https://ropensci.github.io/phruta/reference/sq.retrieve.direct.md)
  where a failed sequence search (caught by `tryCatch`) could crash the
  function with an “argument is of length zero” error instead of just
  skipping that search.
- Fixed a bug in
  [`acc.retrieve()`](https://ropensci.github.io/phruta/reference/acc.retrieve.md)
  where running on an operating system other than macOS, Linux, or
  Windows left the results object unassigned instead of falling back to
  a working code path.
- Removed the duplicated curation logic in
  [`sq.curate()`](https://ropensci.github.io/phruta/reference/sq.curate.md).
  The function now shares one internal implementation between its
  file-based and object-based modes, so a fix only needs to be made once
  going forward.
- Removed a broken `Release` job from the R-CMD-check GitHub Action
  workflow. It referenced a job id that didn’t exist and called scripts
  (`deploy.sh`, a Docker Hub webhook) that aren’t part of this
  repository.
- Added a Bioconductor install step to the R-CMD-check GitHub Action so
  `Biostrings`, `msa`, `odseq`, and `DECIPHER` install correctly on CI.
  These are Bioconductor packages and were never reliably resolved by
  [`remotes::install_deps()`](https://remotes.r-lib.org/reference/install_deps.html)
  alone.
- Updated `actions/checkout`, `actions/cache`, and
  `actions/upload-artifact` in the R-CMD-check GitHub Action to v4.
  GitHub deprecated `actions/cache` v1/v2 and started failing runs that
  use them, and `actions/upload-artifact` was pinned to the floating
  `@main` tag, which can change without warning.
- Applied the same `actions/checkout`/`actions/cache` v4 update to the
  pkgdown and test-coverage GitHub Actions, and added a matching
  Bioconductor install step to both, since both install the package and
  need `Biostrings`, `msa`, `odseq`, and `DECIPHER` to do so.
- Fixed a malformed tag trigger in the pkgdown workflow (`-'*'` instead
  of `- '*'`), which is invalid YAML list syntax and likely meant the
  release-tag trigger was never actually firing.
- Confirmed `Remotes:` in DESCRIPTION still needs `cran/ips`, since
  `ips` was archived on CRAN in November 2025. `Rogue` needs no remote,
  it remains on CRAN.

## `phruta` 0.1.3

- Version accepted in ROpenSci
- `salphycon` released (v. 0.1)
