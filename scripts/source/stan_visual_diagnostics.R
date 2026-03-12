viz_diagnose <- function(fit, pars, warmup=FALSE) {
  if(warmup==FALSE) {
    inc_warmup <- FALSE
    warmup <- 0
  } else inc_warmup <- TRUE

  traceplot <-
    mcmc_trace(
      fit$draws(inc_warmup=inc_warmup),
      pars=grep(pars, names(tidy_draws(fit)), value=TRUE, perl=TRUE),
      n_warmup=warmup
    )
  trankplot <-
    mcmc_rank_overlay(
      fit$draws(),
      pars=grep(pars, names(tidy_draws(fit)), value=TRUE, perl=TRUE)
    )

  traceplot | trankplot
}
