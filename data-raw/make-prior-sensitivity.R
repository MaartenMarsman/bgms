# ==============================================================================
# Generator for vignettes/prior-sensitivity-ps.rds
# ==============================================================================
# The vignette shows a real prior_sensitivity_check() report, not pasted text
# (F-086). Running the check at build time cost the vignette about a minute of
# CPU on every build, on CRAN's machines as well as ours, for an analysis whose
# answer never changes. Following the package's earlier vignettes, the analysis
# runs here, once, locally, and the vignette ships and prints the object.
#
# Run from the package root after any change that could move the report:
#
#     Rscript data-raw/make-prior-sensitivity.R
#
# The script is tracked so the vignette's number has a provenance; dev/ is
# .Rbuildignore'd, so it does not ship. The .rds does ship, and holds the ps
# object and nothing else.
#
# The seeds and arguments below MUST stay identical to the call the vignette
# displays -- that displayed call is the claim that this object came from it.
# ==============================================================================

library(bgms)

data = Wenchuan

fit = bgm(data[, 1:6],
  seed = 1234, chains = 2, cores = 2,
  display_progress = "none", verbose = FALSE
)
ps = prior_sensitivity_check(fit,
  seed = 1234, iter = 3000, warmup = 1000,
  cores = 2
)

# prior_sensitivity_check() keeps the anchor refits only when asked to
# (keep_fits), so at the vignette's call there are no chains attached; the
# assignment below is a belt-and-braces drop, and the size check after the write
# is the rail this actually has to meet.
ps$fits = NULL

out = file.path("vignettes", "prior-sensitivity-ps.rds")
saveRDS(ps, out, version = 2)

size_kb = file.size(out) / 1024
cat(sprintf("wrote %s (%.1f KB)\n", out, size_kb))
if(size_kb >= 100) {
  stop(
    "The vignette payload must stay in double-digit KB; got ",
    round(size_kb, 1), " KB. Thin it further before shipping."
  )
}

# The vignette prints the object; if that errors here it will error at build.
print(ps)
