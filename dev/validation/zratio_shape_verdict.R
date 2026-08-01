# Surface-vs-gold verdict at one Gamma diagonal shape, following the WP4
# scoring protocol: build the absolute-moment surface from real anchors at two
# build seeds, score its deployed per-edge log-ratio against block-Gibbs gold
# over the interior bands, and spot-check the boundary-slope tail past the
# anchored range. The additive kernel it replaces is scored on the same blocks.
#
# Everything runs through the shipped code -- `zratio_build_surfaces` for the
# build and the C++ `zratio_test_surface_eval` for the deploy -- so no part of
# the path being validated is reimplemented here. The only thing this script
# overrides is the deployment fence, which is what the run exists to inform.
#
# Every table carries an already-certified control shape (default 2, scored in
# WP4). Its gold is already banked, so the control costs a surface build and
# reproduces a recorded result or fails the run.
#
# Usage, from the repo root:
#   Rscript dev/validation/zratio_shape_verdict.R --shape=10 --cores=12
#
# Options: --shape (required), --control (default 2, "none" to drop),
#          --cores, --out, --cap, --seeds (build seeds).

suppressMessages(devtools::load_all(".", quiet = TRUE))
source("dev/validation/zratio_gold_bank.R")

opt = function(name, default = NULL) {
  a = grep(paste0("^--", name, "="), commandArgs(TRUE), value = TRUE)
  if(length(a) == 0L) default else sub(paste0("^--", name, "="), "", a[1])
}

shape = as.numeric(opt("shape"))
if(!is.finite(shape)) stop("pass --shape=<value>", call. = FALSE)
control = opt("control", "2")
cores = as.integer(opt("cores", "12"))
cap = as.integer(opt("cap", "80"))
build_seeds = as.integer(strsplit(opt("seeds", "700000,900000"), ",")[[1]])
out_path = opt("out", sprintf("dev/validation/wp5_verdict_%s.rds",
                              sub("[.]", "p", format(shape))))

# The WP4 cell, and its gold budget. Gold seeds and budget match the banked
# {3, 5} rows exactly so a re-score is a lookup and the shapes are comparable.
DELTA = 0.5 * log(12)
ETAS = as.numeric(strsplit(opt("etas", "1,2"), ",")[[1]])
GOLD_SWEEPS = as.integer(opt("sweeps", "4000"))
GOLD_BURN = as.integer(opt("burn", "500"))
# Gold replicate seeds, per shape. The banked references were paid for under
# two vintages -- WP4 scored shapes 0.5 and 2 at seeds 11/23/37, the {3, 5} run
# at seeds 4/8/12 -- so the control shape reuses its own banked seeds and stays
# a lookup instead of buying a third set of replicates for a cell that is
# already certified.
GOLD_SEEDS = as.integer(strsplit(opt("gold_seeds", "4,8,12"), ",")[[1]])
CONTROL_SEEDS = as.integer(strsplit(opt("control_seeds", "11,23,37"), ",")[[1]])

# Interior bands and the far-field spot check. The interior sizes sit well
# inside the anchored range; the far size sits past it, so it scores the
# boundary-slope extension rather than the fit.
SIZES = list(cn = c(20, 26, 32, 38, 42), bip = c(24, 30, 36, 42))
FAR_SIZE = as.numeric(opt("far", "100"))

# A smoke run shrinks the bands and the budget. Its gold keys differ from the
# banked ones, so point BGMS_GOLD_BANK at a scratch file before using it.
if(identical(opt("smoke", "no"), "yes")) {
  SIZES = list(cn = c(20), bip = c(24))
  if(!nzchar(Sys.getenv("BGMS_GOLD_BANK"))) {
    stop("a smoke run must set BGMS_GOLD_BANK to a scratch path", call. = FALSE)
  }
}

control_shape = if(identical(control, "none")) NA_real_ else as.numeric(control)
shapes = unique(c(shape, if(is.finite(control_shape)) control_shape))
seeds_for = function(alpha) {
  if(is.finite(control_shape) && isTRUE(all.equal(alpha, control_shape))) {
    CONTROL_SEEDS
  } else {
    GOLD_SEEDS
  }
}

cat(sprintf(
  "shape verdict: shapes %s | etas %s | cap %d | build seeds %s | cores %d\n",
  paste(shapes, collapse = ", "), paste(ETAS, collapse = ", "), cap,
  paste(build_seeds, collapse = ", "), cores
))
cat(sprintf("gold: %d sweeps + %d burn, seeds %s (control %s: seeds %s)\n\n",
            GOLD_SWEEPS, GOLD_BURN, paste(GOLD_SEEDS, collapse = ", "),
            control, paste(CONTROL_SEEDS, collapse = ", ")))

# The deployment fence gates shapes outside the validated range, which is
# exactly what this run is scoring. Widen it for the duration; the numbers
# decide where it belongs afterwards.
ns = asNamespace("bgms")
fence_old = c(get(".zratio_surface_shape_lo", ns), get(".zratio_surface_shape_hi", ns))
for(nm in c(".zratio_surface_shape_lo", ".zratio_surface_shape_hi")) {
  try(unlockBinding(nm, ns), silent = TRUE)
}
assign(".zratio_surface_shape_lo", 0, envir = ns)
assign(".zratio_surface_shape_hi", Inf, envir = ns)
on.exit({
  assign(".zratio_surface_shape_lo", fence_old[1], envir = ns)
  assign(".zratio_surface_shape_hi", fence_old[2], envir = ns)
}, add = TRUE)

options(bgms.zratio_surface_cache = FALSE, bgms.correction_table_cache = FALSE,
        bgms.verbose = FALSE)

rows = list()
for(alpha in shapes) {
  for(eta in ETAS) {
    # The constants layer is certified separately and first; a shape outside
    # its tuned range warns, which is expected here and not a failure.
    zc = withCallingHandlers(
      bgms:::zratio_constants(DELTA, eta, alpha = alpha),
      warning = function(w) {
        cat("  constants note:", conditionMessage(w), "\n")
        invokeRestart("muffleWarning")
      }
    )
    cell = list(delta = DELTA, eta = eta, alpha = alpha, slab = "normal")

    t0 = proc.time()[["elapsed"]]
    surfs = lapply(build_seeds, function(sd0) {
      bgms:::zratio_build_surfaces(zc, max_size = cap, cores = cores,
                                   seed0 = as.integer(sd0))
    })
    if(any(vapply(surfs, is.null, logical(1)))) {
      stop("surface build returned NULL at shape ", alpha, ", eta ", eta,
           call. = FALSE)
    }
    cat(sprintf("[alpha %g eta %g] %d surface builds in %.0fs\n", alpha, eta,
                length(surfs), proc.time()[["elapsed"]] - t0))

    for(family in c("cn", "bip")) {
      sizes = c(SIZES[[family]], FAR_SIZE)
      gold = gold_bank(zc, cell, family, sizes, GOLD_SWEEPS, GOLD_BURN,
                       seeds_for(alpha), verbose = TRUE, cores = cores)
      for(size in sizes) {
        G = gold_block(family, size)
        gm = gold$mean[gold$size == size]
        gsd = gold$sd[gold$size == size]
        errs = vapply(surfs, function(sf) {
          sv = bgms:::zratio_test_surface_eval(
            G, 1, 2, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
            sf, zc$delta, zc$eta, FALSE, alpha
          )
          abs(sv$logR - gm)
        }, numeric(1))
        add = bgms:::zratio_test_eval(
          G, matrix(c(1, 2), 1, 2), zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt,
          zc$psi0
        )$log_zratio[1]
        rows[[length(rows) + 1L]] = data.frame(
          alpha = alpha, eta = eta, family = family, size = size,
          band = if(size > cap) "far" else "in-hull",
          surface = stats::median(errs), surface_max = max(errs),
          additive = abs(add - gm), gold = gm, gold_sd = gsd
        )
      }
    }
  }
}

verdict = do.call(rbind, rows)
verdict = verdict[order(verdict$alpha, verdict$eta, verdict$family,
                        verdict$size), ]
rownames(verdict) = NULL
saveRDS(verdict, out_path)
gold_bank_manifest()

cat("\n=== verdict ===\n")
print(verdict, digits = 3)
cat("\nwritten to", out_path, "\n")

# Per-cell summary against the recorded 0.003-nat in-range envelope.
ih = verdict[verdict$band == "in-hull", ]
if(nrow(ih) > 0L) {
  cat("\n=== per-cell summary (in-hull band) ===\n")
  agg = stats::aggregate(cbind(surface_max, additive) ~ alpha + eta + family,
                         data = ih, FUN = max)
  agg$inside_envelope = agg$surface_max <= 0.003
  print(agg, digits = 3)

  # Does the reference itself move? A band whose gold is flat in size cannot
  # discriminate between a surface that tracks it and one that returns a
  # constant, so the comparison would have no power. Standing gate in this
  # harness since the WP5 rung-1 reference that never moved.
  cat("\n=== reference movement across the in-hull band ===\n")
  mv = do.call(rbind, lapply(split(ih, list(ih$alpha, ih$eta, ih$family),
                                   drop = TRUE), function(d) {
    span = max(d$gold) - min(d$gold)
    data.frame(
      alpha = d$alpha[1], eta = d$eta[1], family = d$family[1],
      gold_lo = min(d$gold), gold_hi = max(d$gold), span = span,
      max_gold_sd = max(d$gold_sd), span_over_sd = span / max(d$gold_sd),
      moves = span > 3 * max(d$gold_sd)
    )
  }))
  rownames(mv) = NULL
  print(mv, digits = 3)
  if(!all(mv$moves)) {
    cat("\nWARNING: the reference is flat across the band in",
        sum(!mv$moves), "cell(s); those rows do not discriminate.\n")
  }
}
far = verdict[verdict$band == "far", ]
if(nrow(far) > 0L) {
  cat("\n=== far-field (size", FAR_SIZE, ") ===\n")
  print(far[, c("alpha", "eta", "family", "surface", "surface_max", "additive",
                "gold", "gold_sd")], digits = 3)
}
