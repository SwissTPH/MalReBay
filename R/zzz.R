.onAttach <- function(libname, pkgname) {
  # Verify CmdStan is installed and reachable. If not, guide the user once.
  tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(e) {
      packageStartupMessage(
        "MalReBay: CmdStan not found. Run cmdstanr::install_cmdstan() once to install it."
      )
    }
  )
}
