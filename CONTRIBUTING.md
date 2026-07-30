# Contributing to bgms

Thank you for your interest in contributing to bgms. This document
covers the practical steps for building, testing, and submitting
changes.

## Ways to contribute

Contributions come in many forms, and all of them are valued:

- **Report a bug** — open an issue and the bug report template will
  load automatically.
- **Suggest a feature** — open an issue and the feature request template
  will load automatically.
- **Improve documentation** — fix a typo, clarify a vignette, or add an
  example.
- **Write tests** — help us cover edge cases and new functionality.
- **Review a pull request** — a second pair of eyes is always welcome.
- **Contribute code** — fix a bug, implement a feature, or improve
  performance.

For small changes (typos, documentation fixes), open a pull request
directly. For larger work (new features, algorithm changes), open an
issue first so we can discuss the approach before you invest time.

## Getting help

If you get stuck or have questions about contributing, open an issue or
contact the maintainer at m.marsman@uva.nl.

## Getting started

```bash
git clone https://github.com/Bayesian-Graphical-Modelling-Lab/bgms.git
cd bgms
```

Build and check the package:

```r
devtools::document()
devtools::check()
```

## Code style

bgms uses `=` for assignment (not `<-`) and omits the space between
`if`/`for`/`while` and `(`. You do not need to memorize these rules:
running the project styler (step 1 below) applies the formatting, and CI
runs `lintr::lint_package()` on every push and pull request
(`.github/workflows/lint.yaml`) to catch violations.

Before committing, run the following in order:

```r
# 1. Enforce code style
source("inst/styler/bgms_style.R")
styler::style_pkg(style = bgms_style)

# 2. Check for lint issues
lintr::lint_package()

# 3. Regenerate documentation (if you changed any roxygen comments)
roxygen2::roxygenise()
```

The pull request template includes a full checklist of these steps.

## Documentation

- **Exported R functions** get full roxygen2 blocks (`@title`,
  `@description`, `@param`, `@return`, `@examples`, `@seealso`,
  `@export`).
- **Internal R functions** use plain `#` comment blocks, not `#'`
  roxygen.
- **C++ headers** use Doxygen `/** */` blocks on classes and methods.
- **C++ implementations** use inline `//` comments only.

When modifying a function signature, update its documentation in the
same commit.

## Tests

Run the test suite with:

```r
devtools::test()
```

Test files live in `tests/testthat/`. Each `test_that()` description
should be an imperative sentence.

## CI

Pull requests are checked automatically by GitHub Actions:

- `R-CMD-check.yaml` — builds and checks on multiple platforms
- `lint.yaml` — runs `lintr::lint_package()` and fails the PR on lint errors
- `test-coverage.yaml` — reports test coverage

## Commit messages

Use the format: `type: short description`. Common types: `feat`,
`fix`, `docs`, `refactor`, `test`.

## Contributors

All contributors are acknowledged in `inst/CONTRIBUTORS.md`.

## Code of Conduct

This project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to bgms,
you agree to abide by its terms.
