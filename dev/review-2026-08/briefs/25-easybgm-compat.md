# Brief 25 — easybgm compatibility report: what bgms 0.2.0.0 changes downstream (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here.
**REPORT-ONLY for the bgms repo, and you do not modify easybgm at all** — the
deliverable is the report that lets the maintainer (or a later agent, in the
easybgm repo) make the changes. Maintainer request: enumerate every bgms
0.2.0.0 change that reaches easybgm, test the CRAN easybgm against the new
bgms, and verify the compatibility shim still does its work.

## Machine budget (standing rule)

A lead measurement job may be running when you start (long compare fits).
Check `ps aux | grep "[e]xec/R"`: while anything is running, do ONLY the
static work (tasks 1–2 need no fits); hold every fit until the machine is
free, then ~4 threads, one fit at a time. Runtime is not a grading criterion.

## Setup

- Repo (Dropbox; do NOT switch its branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Base on **`origin/develop`** (NOT the local ref) AT OR AFTER `a37c045d`:
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-val25 -b review/easybgm-compat origin/develop
  cd ~/bgms-review/wt-val25
  ```
- Private library `~/bgms-review/lib-val25`: install this bgms build into it,
  plus **easybgm from CRAN** (binary is fine; record the exact version) and
  its dependencies. Also download the easybgm **source tarball** from CRAN —
  task 2 reads its code, not its docs.
- Scratch under `~/bgms-review/val25/`. The ONLY tracked changes on your
  branch are the report. No pushing; no attribution trailers. Everything in
  `R/`, `src/`, `tests/`, `man/`, `vignettes/`, `NEWS.md` is off-limits.

## Tasks

### 1. The compatibility shim, mapped then exercised

bgms carries a version gate for easybgm: `needs_easybgm_s3_compat()`
(`R/build_output.R:307-326`) and a deprecated easybgm-compatibility section
in the S7 class (`R/class_s7.R:92` onward). Map what the gate does in each of
the three states — easybgm not loaded / easybgm < 0.5.0 loaded / easybgm
>= 0.5.0 loaded — with file:line citations, including exactly what object
shape each state hands easybgm. Then EXERCISE it live against the CRAN
easybgm: load it, fit a small `bgm()`, and show which branch fires and that
its behaviour is the designed one (clean extractor path, or the instructive
error verbatim). If the CRAN version sits below 0.5.0, that error IS the
designed behaviour — demonstrate it and say so plainly.

### 2. The delta enumeration — every 0.2.0.0 change that reaches easybgm

`NEWS.md` is the authoritative delta record against 0.1.6.3 (the version
CRAN easybgm was built against). For every entry that could reach easybgm,
read easybgm's CRAN source and answer whether it actually does — cite
easybgm file:line for every claim. Classify each item:
**BREAKS** (easybgm errors or misbehaves) / **CHANGES NUMBERS** (same code
runs, different results a user or test would notice) / **COSMETIC** /
**NOT REACHED**. Cover at least:

- the S7 fit object and the extractor contract (which easybgm accessors
  survive; what the class_s7 compat section still provides and that it is
  deprecated);
- defaults that move numbers: `bgm()` interaction prior Cauchy(2.5-era) →
  `normal_prior(1)`; `bgmCompare()` baseline AND `difference_family` → Normal;
  `precision_scale_prior`'s eta frame; `precision_graph_prior = "hierarchical"`
  default (GGM/mixed only); `iter`/`warmup` = 2e3;
- removed/deprecated formals easybgm might pass (`standardize`,
  `pairwise_scale`, `interaction_scale`, `main_alpha`/`main_beta`, ...);
- `extract_ess()` indicator semantics (RB n_eff now), `extract_inclusion_bf`'s
  `log=`, natural-log BF display conventions;
- compare-side union category semantics (`category_support`, the classed
  group-support warning, the `*` marking in summaries) if easybgm touches
  compare output at all;
- anything else your NEWS read turns up that easybgm's source reaches.

### 3. Run the CRAN easybgm against this bgms (fits — machine-free only)

easybgm's own check machinery first: its test suite run against the new bgms
build (note which tests even reach bgms). Then one representative
easybgm-workflow per model family easybgm exposes through bgms (small
`iter`), reading its vignette/examples for the canonical calls. Record every
error, warning, and changed-output class verbatim. This is the dress
rehearsal for Phase 3's reverse-dependency check — report what
`R CMD check` of easybgm against the new bgms would show, and run exactly
that if the machine allows it (easybgm's check is small).

### 4. The recommendations — the part the maintainer asked for

- Per-file, per-line recommended changes in easybgm (or "no change needed —
  numbers shift only"), written so an agent in the easybgm repo could execute
  them without reading bgms first.
- A short draft notification memo (the downstream-maintainer message: what
  changed, what to do, in ten lines) — the lead folds it into the Phase-3
  memo.
- **The ordering question, answered explicitly:** does bgms 0.2.0.0 break the
  CRAN easybgm outright? If yes, CRAN's reverse-dependency check will flag
  it, and the report must say what the submission order should be
  (easybgm update first / simultaneous / bgms first with notification), with
  the evidence from task 3.

## Verification gate

1. Every claim about easybgm cites its source (file:line in the CRAN
   tarball); every BREAKS/CHANGES-NUMBERS classification is backed by an
   executed demonstration where one is possible.
2. The shim's live branch demonstration (task 1) is in the report verbatim.
3. Task-3 outcomes reported for every workflow attempted; nothing dropped
   silently.
4. easybgm's exact CRAN version and the bgms commit recorded.

## Deliverable

`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/25-easybgm-compat.md`:
Shim map + live demonstration / the classified delta table with citations /
task-3 run log summary / per-file recommendations + draft memo + the ordering
answer / findings (severity-tagged) / open questions. Copy to the Dropbox
path; commit the report on the branch.
