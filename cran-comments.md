# bgms 0.2.0.0 — CRAN submission comments

## Test environments

- Local: macOS 15 (aarch64), R 4.6.0 — `R CMD check --as-cran`
- win-builder (release, devel) — pending
- mac-builder — pending

## R CMD check results

0 ERRORs, 0 WARNINGs. Local NOTEs are environment artifacts only (HTML Tidy
version; incoming-feasibility remarks for a version bump).

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
