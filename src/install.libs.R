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
#
# PATCHED (2026-09-10): skip recompilation for any model that already
# ships a precompiled executable in src/stan/ (checked into the repo).
# CmdStan/a C++ toolchain aren't available on managed deploy targets like
# shinyapps.io, so those models must already be built and portable
# (relative $ORIGIN rpath, bundled .so deps) before install ever runs.
# Models with no shipped executable still compile normally here, so local
# development on a machine with CmdStan installed is unaffected.
#
# PATCHED (2026-10-04): the shipped executable is a Linux x86-64 binary, so
# it is only used on Linux x86-64 (e.g. shinyapps.io); macOS, Windows and
# other architectures compile from source. It is also only used when the
# .stan file still matches the MD5 fingerprint in <model>.stan.md5, written
# when the binary was built -- so editing the model without rebuilding the
# binary recompiles instead of silently shipping the old model. After
# rebuilding the binary, update the fingerprint (from the package root):
#   writeLines(unname(tools::md5sum("src/stan/malrebay_model.stan")),
#              "src/stan/malrebay_model.stan.md5")
# Finally, the binary must actually run here (`<model> info` loads its
# shared libraries and exits in milliseconds): on a system whose glibc /
# libstdc++ are older than the build machine's, it recompiles -- or, with no
# CmdStan (shinyapps.io), fails the install now rather than the first
# analysis in the deployed app.
# -----------------------------------------------------------------------
callr::r(
  func = function(bin_stan) {
    linux_x86_64 <- Sys.info()[["sysname"]] == "Linux" &&
      Sys.info()[["machine"]] %in% c("x86_64", "amd64")
    runs_here <- function(exe) {
      status <- tryCatch(
        suppressWarnings(system2(exe, "info", stdout = FALSE, stderr = FALSE)),
        error = function(e) 1L
      )
      identical(as.integer(status), 0L)
    }
    models <- instantiate::stan_package_model_files(path = bin_stan)
    to_compile <- Filter(function(model) {
      exe <- tools::file_path_sans_ext(model)
      if (.Platform$OS.type == "windows") exe <- paste0(exe, ".exe")
      if (!file.exists(exe)) return(TRUE)
      fingerprint <- paste0(model, ".md5")
      matches_source <- file.exists(fingerprint) &&
        identical(trimws(readLines(fingerprint, n = 1, warn = FALSE)),
                  unname(tools::md5sum(model)))
      usable <- linux_x86_64 && matches_source &&
        file.access(exe, mode = 1) == 0
      precompiled <- usable && runs_here(exe)
      if (precompiled) {
        message("Using precompiled Stan model: ", exe)
      } else {
        if (linux_x86_64 && !matches_source) {
          message("Precompiled Stan model is out of date with ", basename(model),
                  " (or has no .md5 fingerprint); recompiling.")
        } else if (usable) {
          message("Precompiled Stan model can't run on this system (missing or ",
                  "too-old system libraries); recompiling.")
        }
        # Remove the unusable binary so cmdstanr can't treat it as up to date
        unlink(exe)
      }
      !precompiled
    }, models)
    for (model in to_compile) {
      tryCatch(
        cmdstanr::cmdstan_model(stan_file = model, compile = TRUE, force_recompile = TRUE),
        error = function(e) stop(
          "Compiling the MalReBay Stan model failed. Installing MalReBay here ",
          "needs CmdStan (see cmdstanr::install_cmdstan()) because the shipped ",
          "precompiled model can't be used on this system.\n",
          "Original error: ", conditionMessage(e), call. = FALSE
        )
      )
    }
  },
  args = list(bin_stan = bin_stan),
  show = TRUE,
  stderr = "2>&1"
)
