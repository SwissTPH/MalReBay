.onAttach <- function(libname, pkgname) {
  # Verify CmdStan is installed and reachable. If not, guide the user once.
  ok <- tryCatch({
    path <- cmdstanr::cmdstan_path()
    nzchar(path) && file.exists(path)
  }, error = function(e) FALSE)
  
  if (!ok) {
    packageStartupMessage(
      "NOTE: CmdStan is not installed or not found.\n",
      "MalReBay requires CmdStan to run. Install it with:\n",
      "  cmdstanr::install_cmdstan()\n",
      "See ?cmdstanr::install_cmdstan for details."
    )
  }
}
