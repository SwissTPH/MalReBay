# Check MCMC sampler diagnostics and report plain-language guidance

After Stan sampling completes for a site, inspects the raw HMC/NUTS
sampler diagnostics (divergent transitions, max-treedepth hits, E-BFMI)
and, for any that look problematic, prints a plain-language explanation
of what it means and what to change in `mcmc_config` to fix it. This is
a supplement to cmdstanr's own diagnostic messages (already printed
automatically during `$sample()`), not a replacement – those are more
technical but this translates them into a concrete next step.

## Usage

``` r
check_mcmc_diagnostics(fit, site_name, adapt_delta, verbose = TRUE)
```

## Arguments

- fit:

  A `CmdStanMCMC` fit object, as returned by `$sample()`.

- site_name:

  A character string used to label messages.

- adapt_delta:

  The `adapt_delta` value used for this fit (echoed back in the
  divergence message so the suggested next value is relative to what was
  actually tried).

- verbose:

  Logical. If `FALSE`, no messages are printed and the diagnostics are
  only returned invisibly.

## Value

Invisibly, the list returned by `fit$diagnostic_summary()`
(`num_divergent`, `num_max_treedepth`, `ebfmi`; one value per chain), or
`NULL` if diagnostics could not be computed.

## Details

- **Divergent transitions** mean the sampler lost numerical accuracy
  while exploring some region of the posterior; results for that site
  may be biased until this is resolved. The fix is to increase
  `adapt_delta` in `mcmc_config` (e.g. to `0.95` or `0.99`), which
  forces smaller, more careful sampling steps at the cost of speed.

- **Hitting the maximum tree depth** is an efficiency issue, not a bias
  one – sampling was simply less thorough in those iterations. On its
  own (no divergences, good R-hat/ESS) it can usually be tolerated;
  alongside divergences, fix `adapt_delta` first.

- **Low or undefined E-BFMI** means a chain explored the tails of the
  posterior poorly, which can make the resulting probabilities too
  overconfident. The fix is to increase `iter` and/or `n_chains` in
  `mcmc_config` and rerun.
