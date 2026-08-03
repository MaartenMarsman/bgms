# bgms 0.2.0.0 — CRAN submission comments

## Test environments

- Local: macOS 15 (aarch64), R 4.6.0 — `R CMD check --as-cran`

## R CMD check results

0 ERRORs, 0 WARNINGs. Local NOTEs are environment artifacts only (HTML Tidy
version; incoming-feasibility remarks for a version bump).

## Notes from the previous submission

- The "(possibly) invalid file URIs" NOTE (CONTRIBUTING.md and
  CODE_OF_CONDUCT.md linked from README.md) is fixed: README.md now links to
  the files on GitHub with absolute URLs, since both files are excluded from
  the tarball via .Rbuildignore.
- The "Overall checktime 13 min > 10 min" NOTE on r-devel-windows-x86_64:
  the four heaviest test files (internal numerical regression guards, ~5 of
  the 8.4 test minutes on that machine) are now skipped on CRAN and run on
  the package's own CI instead (every push, plus nightly and weekly tiers).
  The product-surface tests still run on CRAN in full.
- The "new files in some other directories" NOTE on the
  r-devel-linux-x86_64-debian-gcc check of the published version (0.1.6.3)
  lists Rtmp/xvfb-run/pocl/shm files from concurrent checks on the shared
  check machine; they are not created by this package. The package's own
  checks report "checking for detritus in the temp directory ... OK", and the
  package creates no files under ~/tmp/scratch, /dev/shm, or ~/.cache.

## Reverse dependencies

easybgm (0.4.0) was checked against this version and against bgms 0.1.6.3:
identical check status (no new ERROR, WARNING, or NOTE), and its test suite
passes 0-failure/0-error against both. A compatibility layer in bgms serves
easybgm <= 0.4.0 the object layout it expects. The easybgm maintainers have
been notified with a per-line change list for their next release; the only
visible change on their side is runtime (bgms's default `warmup` doubled).

## Notes for the reviewers

- This release replaces the samplers wholesale. On seeded synthetic datasets
  where 0.1.6.3's ordinal sampler fails to converge (split-R-hat up to 12.3
  and effective sample sizes near zero on three of five seeds), this version
  converges on all five and lands closer to the generating truth on all five.
- The statistical certification suites (simulation-based calibration,
  parameter recovery, cross-method validation) are tiered behind environment
  variables and run on scheduled CI rather than inside `R CMD check`; the
  check runs the fast deterministic tier (~8,500 expectations, ~90 s).
- Examples and vignettes are seeded and runtime-capped for CRAN machines.
