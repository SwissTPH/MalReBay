libs <- file.path(R_PACKAGE_DIR, "libs", R_ARCH)
dir.create(libs, recursive = TRUE, showWarnings = FALSE)
for (file in c("symbols.rds", Sys.glob(paste0("*", SHLIB_EXT)))) {
  if (file.exists(file)) {
    file.copy(file, file.path(libs, file))
  }
}
inst_stan <- file.path("..", "inst", "stan")
if (dir.exists(inst_stan)) {
  warning(
    "Stan models in inst/stan/ are deprecated in {instantiate} ",
    ">= 0.0.4.9001 (2024-01-03). Please put them in src/stan/ instead."
  )
  if (file.exists("stan")) {
    warning("src/stan/ already exists. Not copying models from inst/stan/.")
  } else {
    message("Copying inst/stan/ to src/stan/.")
    fs::dir_copy(path = inst_stan, new_path = "stan")
  }
}
bin <- file.path(R_PACKAGE_DIR, "bin")
if (!file.exists(bin)) {
  dir.create(bin, recursive = TRUE, showWarnings = FALSE)
}
bin_stan <- file.path(bin, "stan")
fs::dir_copy(path = "stan", new_path = bin_stan)
# -----------------------------------------------------------------------
# PATCHED: compile directly via cmdstanr instead of
# instantiate::stan_package_compile(). That helper currently hardcodes a
# `threads=` argument when calling cmdstanr's $compile() method, which
# recent cmdstanr versions (moving toward v1.0) no longer accept, failing
# with "unused argument (threads = FALSE)".
# See https://github.com/wlandau/instantiate/issues/37 (unresolved as of
# 2026-07-13). This block only uses cmdstanr::cmdstan_model(stan_file=,
# compile=TRUE), which is stable across cmdstanr versions and produces the
# executable at the same path instantiate::stan_package_model() expects to
# find at runtime (same directory as the .stan file, same base name, plus
# ".exe" on Windows — cmdstanr's own default naming, so no exe_file= override
# is needed here).
# TODO: once instantiate fixes #37, revert to:
#   instantiate::stan_package_compile(
#     models = instantiate::stan_package_model_files(path = bin_stan)
#   )
# -----------------------------------------------------------------------
callr::r(
  func = function(bin_stan) {
    models <- instantiate::stan_package_model_files(path = bin_stan)
    for (model in models) {
      cmdstanr::cmdstan_model(stan_file = model, compile = TRUE)
    }
  },
  args = list(bin_stan = bin_stan),
  show = TRUE,
  stderr = "2>&1"
)
