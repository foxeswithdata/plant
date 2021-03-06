##' Build an appropriately refined schedule.
##'
##' There are control options (within the \code{Parameters} object)
##' that affect how this function runs, in particular
##' \code{schedule_nsteps} and \code{schedule_eps} control how refined
##' the schedule will end up, and \code{schedule_verbose} controls if
##' details are printed to the screen during construction.
##'
##' @title Build Cohort Schedule
##' @param p Parameters object
##' @return A Parameters object, with schedule components set.  The
##' output seed rain is also available as an attribute
##' \code{seed_rain}.
##' @author Rich FitzJohn
##' @export
build_schedule <- function(p) {
  p <- validate(p)

  n_spp <- length(p$strategies)
  if (n_spp == 0L || !any(p$is_resident)) {
    stop("Can't build a schedule with no residents")
  }
  control <- p$control
  eps <- control$schedule_eps

  for (i in seq_len(control$schedule_nsteps)) {
    res <- run_scm_error(p)
    seed_rain_out <- res[["seed_rain"]]
    split <- lapply(res$err$total, function(x) x > eps)

    if (!any(unlist(split), na.rm=TRUE)) {
      break
    }

    ## Prepare for the next iteration:
    times <- p$cohort_schedule_times
    for (idx in seq_len(n_spp)) {
      times[[idx]] <- split_times(times[[idx]], split[[idx]])
    }
    p$cohort_schedule_times <- times

    msg <- sprintf("%d: Splitting {%s} times (%s)",
                   i,
                   paste(sapply(split, sum),    collapse=","),
                   paste(sapply(split, length), collapse=","))
    plant_log_debug(msg, routine="schedule", event="split", round=i)
  }

  p$cohort_schedule_ode_times <- res$ode_times
  ## Useful to record the last seed rain out:
  attr(p, "seed_rain_out") <- seed_rain_out

  p
}

split_times <- function(times, i) {
  ## Upwind splitting scheme only, which means that we will never
  ## split the last interval [assuming OK for now].  Inefficiently
  ## interleaved with sort().  These issues can change easily enough
  ## later.  The aim is making sure that we don't introduce the same
  ## point twice; one from upstream and one from downstream.
  dt <- diff(times)
  i <- which(i)
  sort(c(times, times[i] - dt[i-1]/2))
}
