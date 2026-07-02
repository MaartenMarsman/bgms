# Quiet bgms's advisory verbose output (warmup-length warnings, data-cleaning
# messages) during the test run. Tests that specifically assert this output set
# options(bgms.verbose = TRUE) locally, which overrides this default.
options(bgms.verbose = FALSE)

# Keep edge-prior correction tables out of the user-level cache directory:
# any test fit with a hierarchical edge prior on continuous or mixed data
# builds one at fit time. Identity tests that share a cell set their own
# directory locally.
options(
  bgms.correction_cache_dir = file.path(tempdir(), "bgms-ctable-tests")
)
