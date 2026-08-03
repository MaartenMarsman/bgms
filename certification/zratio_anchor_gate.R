# WP4 attribution gate and cap-80 build cost, re-derived under the slice-pivot
# anchor kernel.
#
# Two questions, in this order:
#
# 1. Anchor sweep multipliers. The gate matches the anchor Monte-Carlo error at
#    a non-unit shape to the shape-1 reference, measured as the across-seed
#    spread of log S1 on a fixed block, rather than inferred from acceptance.
#    Under the independence-Metropolis kernel this forced 2x (4x at eta 2,
#    cn, shape 0.5); the slice kernel's pivot ESS is flat in shape and size, so
#    the multipliers are re-derived here rather than carried over.
#
# 2. Cap-80 build cost, serial, against the standing <= 60 s acceptance.
#
# TIMINGS REQUIRE AN INSTALLED BUILD. devtools::load_all compiles at -O0, so a
# build cost measured under it is not a benchmark. This script therefore uses
# the installed package; run `R CMD INSTALL .` first and do not run it
# alongside other load.
#
# Usage, from the repo root:
#   Rscript certification/zratio_anchor_gate.R --cores=12
#   Rscript certification/zratio_anchor_gate.R --cost-only --cores=1

suppressMessages(library(bgms))

opt = function(name, default = NULL) {
  a = grep(paste0("^--", name, "="), commandArgs(TRUE), value = TRUE)
  if(length(a) == 0L) default else sub(paste0("^--", name, "="), "", a[1])
}
has_flag = function(name) any(commandArgs(TRUE) == paste0("--", name))

cores = as.integer(opt("cores", "12"))
DELTA = 0.5 * log(12)
ETAS = c(1, 2)
SHAPES = as.numeric(strsplit(opt("shapes", "0.5,2,3,5,10"), ",")[[1]])
SEEDS = as.integer(opt("seeds", "8"))
GATE_N = 20L
GATE_DENS = 0.9
GATE_SWEEPS = 1000L
MULT_LADDER = c(1L, 2L, 4L)   # hard cap: no search past 4x

ns = asNamespace("bgms")
fence_old = c(get(".zratio_surface_shape_lo", ns),
              get(".zratio_surface_shape_hi", ns))
for(nm in c(".zratio_surface_shape_lo", ".zratio_surface_shape_hi")) {
  try(unlockBinding(nm, ns), silent = TRUE)
}
assign(".zratio_surface_shape_lo", 0, envir = ns)
assign(".zratio_surface_shape_hi", Inf, envir = ns)

options(bgms.zratio_surface_cache = FALSE, bgms.correction_table_cache = FALSE,
        bgms.verbose = FALSE)

constants = function(eta, alpha) {
  suppressWarnings(bgms:::zratio_constants(DELTA, eta, alpha = alpha))
}

# Across-seed spread of log S1 on one fixed block at one sweep budget. This is
# the anchor Monte-Carlo error the surface fit inherits.
anchor_spread = function(eta, alpha, family, sweeps, cores) {
  zc = constants(eta, alpha)
  fn = if(family == "cn") bgms:::zratio_anchor_cn else bgms:::zratio_anchor_bip
  run = function(sd_) {
    a = fn(GATE_N, GATE_DENS, zc, as.integer(sweeps), 200L, as.integer(sd_))
    if(is.null(a)) NA_real_ else log(a$S1)
  }
  seeds = 5000L + seq_len(SEEDS)
  vals = if(cores > 1L && .Platform$OS.type == "unix") {
    unlist(parallel::mclapply(seeds, run, mc.cores = cores,
                              mc.preschedule = FALSE))
  } else {
    vapply(seeds, run, numeric(1))
  }
  stats::sd(vals)
}

if(!has_flag("cost-only")) {
  cat("=== attribution gate: anchor spread vs the shape-1 reference ===\n")
  cat(sprintf("block n=%d density=%g, %d sweeps nominal, %d seeds\n\n",
              GATE_N, GATE_DENS, GATE_SWEEPS, SEEDS))
  rows = list()
  for(eta in ETAS) {
    for(family in c("cn", "bip")) {
      ref = anchor_spread(eta, 1, family, GATE_SWEEPS, cores)
      for(alpha in SHAPES) {
        # Walk the ladder until the spread matches the shape-1 reference.
        # Parity is a ratio at or below 1.2: the reference is itself an
        # 8-seed estimate, so demanding equality would chase its own noise.
        chosen = NA_integer_
        detail = list()
        for(m in MULT_LADDER) {
          s = anchor_spread(eta, alpha, family, GATE_SWEEPS * m, cores)
          detail[[as.character(m)]] = s / ref
          if(s / ref <= 1.2) {
            chosen = m
            break
          }
        }
        rows[[length(rows) + 1L]] = data.frame(
          eta = eta, family = family, alpha = alpha, ref_sd = ref,
          ratio_1x = detail[["1"]],
          ratio_2x = if(is.null(detail[["2"]])) NA_real_ else detail[["2"]],
          ratio_4x = if(is.null(detail[["4"]])) NA_real_ else detail[["4"]],
          resolved = chosen,
          shipped = bgms:::zratio_anchor_shape_multiplier(alpha)
        )
        cat(sprintf("  eta %g %s shape %-4g ref_sd %.4g  1x ratio %.2f -> %s\n",
                    eta, family, alpha, ref, detail[["1"]],
                    if(is.na(chosen)) "NOT MET at 4x" else paste0(chosen, "x")))
        utils::flush.console()
      }
    }
  }
  gate = do.call(rbind, rows)
  saveRDS(gate, "certification/wp5_anchor_gate.rds")
  cat("\n=== resolved multipliers ===\n")
  print(gate, digits = 3)
  cat("\nwritten to certification/wp5_anchor_gate.rds\n")
}

if(has_flag("gate-only")) {
  assign(".zratio_surface_shape_lo", fence_old[1], envir = ns)
  assign(".zratio_surface_shape_hi", fence_old[2], envir = ns)
  quit(save = "no")
}

# The cost figures are benchmarks, so they need a quiet machine and the
# installed (-O2) build; run this arm on its own.
cat("\n=== cap-80 build cost, serial (installed build) ===\n")
cost_shapes = as.numeric(strsplit(opt("cost_shapes", "1,10"), ",")[[1]])
crows = list()
for(alpha in cost_shapes) {
  zc = constants(2, alpha)
  t0 = proc.time()[["elapsed"]]
  surf = bgms:::zratio_build_surfaces(zc, max_size = 80L, cores = 1L)
  el = proc.time()[["elapsed"]] - t0
  crows[[length(crows) + 1L]] = data.frame(
    alpha = alpha, eta = 2, cap = 80L, serial_s = el,
    mult = bgms:::zratio_anchor_shape_multiplier(alpha),
    ok = !is.null(surf), within_60s = el <= 60
  )
  cat(sprintf("  shape %-4g multiplier %dx -> %.1f s serial %s\n", alpha,
              bgms:::zratio_anchor_shape_multiplier(alpha), el,
              if(el <= 60) "(inside 60 s)" else "(OVER 60 s)"))
  utils::flush.console()
}
cost = do.call(rbind, crows)
saveRDS(cost, "certification/wp5_build_cost.rds")
print(cost, digits = 3)

assign(".zratio_surface_shape_lo", fence_old[1], envir = ns)
assign(".zratio_surface_shape_hi", fence_old[2], envir = ns)
