#' Multi-epoch Monte Carlo SEIR
#'
#' Conducts Monte Carlo simulations for a multi-epoch SEIR model
#'
#' @param n.rep The number of simulation replicates. Defaults to 10.
#' @param paramsList A list of parameter value lists. Each list entry must match \code{params} from \code{mcSeirScenario}.
#' @param init Initial values for \code{S}, \code{E}, \code{I}, \code{R} and if \code{func = seira}, \code{A}, the ambient source of infection.
#' @param timesList A list of time vectors.  Each list entry must match \code{times} from \code{mcSeirScenario} and be consecutive.
#' @param interventionList A list of intervention vectors. Each list entry must match \code{intervention} from \code{mcSeirScenario}.
#' @param nThreads The number of threads to use in multi-core computing.
#' @param func The function defining the derivatives in the ODE system (which define the model). Defaults to \code{seira}.
#' @importFrom tidyselect all_of
#' @export
mcMultiEpochSeirScenario <- function(
  n.rep = 10,
  paramsList,
  timesList,
  interventionList,
  init,
  nThreads = 1,
  func = seira,
  modelOutput = if(is.list(init) & !is.data.frame(init)) names(init[[1]]) else names(init)
) {
  if (length(paramsList) != length(timesList) ||
      length(timesList) != length(interventionList)) {
    stop("Must specify the same number of epoches in `paramsList`, `timesList` and `interventionList`")
  }

  transformedTimesList <- transformEpochTimes(timesList)

  compartmentNames <- names(init)

  # First epoch
  result <- mcSeirScenario(n.rep, paramsList[[1]], init, transformedTimesList[[1]],
                           interventionList[[1]], nThreads, func, model.output = modelOutput)

  last <- result
  if (length(paramsList) > 1) {
    for (e in 2:length(paramsList)) {
      df <- last$out %>% filter(.data$times == max(.data$times)) %>% arrange(replicate) %>% select(tidyselect::all_of(compartmentNames))
      init <- split(df, seq(nrow(df)))
      epoch <- mcSeirScenario(n.rep, paramsList[[e]], init, transformedTimesList[[e]],
                              interventionList[[e]], nThreads, func, model.output = modelOutput)
      epoch$out$times <- epoch$out$times + min(timesList[[e]])
      last <- epoch
      result$out <- rbind(result$out, epoch$out[-1,])
    }
  }

  return(result)
}

transformEpochTimes <- function(timesList) {

  transformedTimeList <- list()
  startTime <- 0
  for (epoch in 1:length(timesList)) {
    times <- timesList[[epoch]]
    if (times[1] != startTime) {
      stop("Must specify consecutive times in `timesList`")
    }
    startTime <- times[length(times)]
    transformedTimeList[[epoch]] <- times - times[1]
  }

  return(transformedTimeList)
}

mcSeirScenario <- function(
  n.rep = 10,
  params,
  init,
  times = seq(0, 365, by = 1),
  intervention = c(1, 1, 1),
  nThreads = 1,
  func = seira,
  model.output = if(is.list(init) & !is.data.frame(init)) names(init[[1]]) else names(init),
  c.n = NULL, c.pw.intv = NULL, c.w = NULL, c.pop = NULL, R0 = NULL, i = NULL, L = NULL
)
{
  n.t <- length(times)
  n.out <- n.t * n.rep

  required.params <- c("R0.mean", "R0.sd", "c.w.mean", "c.w.sd",
                       "c.n.mean", "c.n.sd", "c.pw.mean", "c.pw.sd",
                       "c.pop.mean", "c.pop.sd", "i.mean", "i.sd",
                       "L.mean", "L.sd", "a", "b", "c")

  for(j in 1:length(required.params)) {
    if(!(required.params[j] %in% names(params))) {
      if(grepl("mean", required.params[j])) {
        params[required.params[j]] <- 1
      } else {
        params[required.params[j]] <- 0
      }
    }
  }

  if (is.null(c.n)) {
    c.n       <- rgammaMeanSd(n.rep, params[["c.n.mean"]],  params[["c.n.sd"]])
  }
  if (is.null(c.pw.intv)) {
    c.pw.intv <- rgammaMeanSd(n.rep, params[["c.pw.mean"]], params[["c.pw.sd"]]) * intervention[3]
  }
  if (is.null(c.w)) {
    c.w       <- rgammaMeanSd(n.rep, params[["c.w.mean"]], params[["c.w.sd"]])
  }
  if (is.null(c.pop)) {
    c.pop     <- rgammaMeanSd(n.rep, params[["c.pop.mean"]], params[["c.pop.sd"]])
  }
  if (is.null(R0)) {
    R0        <- rgammaMeanSd(n.rep, params[["R0.mean"]], params[["R0.sd"]])
  }
  if (is.null(i)) {
    i         <- rgammaMeanSd(n.rep, params[["i.mean"]], params[["i.sd"]])
  }
  if (is.null(L)) {
    L         <- rgammaMeanSd(n.rep, params[["L.mean"]], params[["L.sd"]])
  }

  R0.intv   <- R0 * intervention[1] * intervention[2] * (c.w/c.pop) # Within-office R0 under intervention
  R0.ext.pop.intv <- R0 * ((c.n + c.pw.intv) / (c.pop)) # office-non office R0 under intervention
  beta  <- R0.intv / i           # beta parameter for employee-employee contact
  gamma <- (1/i)                 # gamma = recovery rate (1/i)
  rho   <- (R0.ext.pop.intv/i)   # "beta" parameter for employee-non-employee contact (set = 0 for traditional SEIR model)

  # Set up output data frame -----------------------------------------------------
  if (is.list(init)) {
    template_init <- init[[1]]
  } else {
    template_init <- init
  }

  seir.out <- data.frame(matrix(NA, ncol = 2 + length(model.output), nrow = n.out))

  if(all.equal(func, seira) == TRUE) {

    if (nThreads == 1) {
      for(j in 1:n.rep) {
        # Writing onto row [(j-1)*n.t+1] to row [j*n.t]:
        seir.out[((j-1)*n.t + (1:n.t)),] <- execute.seira(j, n.t, beta, gamma, rho, L,
                                                          params[["a"]], params[["b"]], params[["c"]],
                                                          init, times, func)[, c("replicate", "time",   model.output)]
      }
    } else {
      cluster <- ParallelLogger::makeCluster(numberOfThreads = nThreads)
      resultList <- ParallelLogger::clusterApply(cluster, 1:n.rep, execute.seira, n.t,
                                                 beta, gamma, rho, L,
                                                 params[["a"]], params[["b"]], params[["c"]],
                                                 init, times, func)
      ParallelLogger::stopCluster(cluster)

      seir.out <- do.call("rbind", resultList)[, c("replicate", "time", model.output)]
    }
  } else if((all.equal(func, sir7) == TRUE)) {

    if (nThreads == 1) {
      for(j in 1:n.rep) {
        # Writing onto row [(j-1)*n.t+1] to row [j*n.t]:
        seir.out[((j-1)*n.t + (1:n.t)),] <- execute.seirmd(j, n.t, params,
                                                           init, times, func)[, c("replicate", "time",   model.output)]
      }
    } else {
      cluster <- ParallelLogger::makeCluster(numberOfThreads = nThreads)
      resultList <- ParallelLogger::clusterApply(cluster, 1:n.rep, execute.seirmd, n.t,
                                                 params, init, times, func)
      ParallelLogger::stopCluster(cluster)

      seir.out <- do.call("rbind", resultList)[, c("replicate", "time", model.output)]
    }

  }
  names(seir.out) <- c("replicate", "times", model.output)

  seir.out$intervention <- attr(intervention, "label")

  # The final output:
  return(list(out = seir.out, R0.internal = data.frame(R0.internal = R0.intv),
              R0.external = data.frame(R0.external = R0.ext.pop.intv),
              R0 = data.frame(R0 = R0.intv + R0.ext.pop.intv),
              i = i, L = L, c.n = c.n, c.w = c.w,
              c.pw.intv = c.pw.intv, R0 = R0,
              label = attr(intervention, "label")))
}

# for each scenario replicate --------------------------------------------------
execute.seira <- function(j, n.t,
                          beta, gamma, rho, L, a, b, c,
                          in_init, times, func) {

  if (is.list(in_init)) {
    init <- in_init[[j]]
    names <- names(init)
    init <- as.numeric(init)
    names(init) <- names
  } else {
    init <- in_init
  }
  parameters <- c(beta = beta[j],
                  gamma = gamma[j],
                  rho = rho[j],
                  L = L[j],
                  a = a,
                  b = b,
                  c = c)     # this is the set of parameters that enters into the differential equations

  # Rate-limiting step
  result <- cbind(
    replicate = rep(j, n.t), # The first column indicates the replicate, which runs from 1 to n.rep
    as.data.frame(deSolve::ode(y = init, times = times, func = func, parms = parameters))
  )

  return(result)
}

execute.seirmd <- function(j, n.t, parameters, in_init, times, func) {

  if (is.list(in_init)) {
    init <- in_init[[j]]
    names <- names(init)
    init <- as.numeric(init)
    names(init) <- names
  } else {
    init <- in_init
  }

  result <- cbind(
    replicate = rep(j, n.t),
    as.data.frame(deSolve::ode(y = init, times = times, func = func, parms = parameters))
  )

  return(result)
}

#' Gamma RNG
#'
#' Generates random samples from a gamma distribution with parameters \code{mean} and \code{sd}.
#'
#' @param n The number of samples.
#' @param mean The mean of the gamma distribution.
#' @param sd The standard deviation of the gamma distribution.
#' @export
rgammaMeanSd <- function(n, mean, sd) {
  if(sd == 0) {
    return(rep(mean, n))
  }
  else {
    return(rgamma(n, shape = mean^2/sd^2, rate = mean/sd^2))
  }
}

#' SEIR Model with Ambient Source of Infections
#'
#' This function defines the system of ODEs for a SEIR model with an ambient source of infections. For use with deSolve::ode.
#'
#' @param times The times at which to evaluate solutions.
#' @param state The current estimate of the variables in the ODE system.
#' @param parameters A vector of parameters. Must include \code{beta}, \code{rho}, \code{L}, \code{gamma}, \code{a}, \code{b} and \code{c}.
#' @export
seira <- function(times, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- - beta * S * I - A * rho * S
    dE <-   beta * S * I + A * rho * S - (E/L)
    dI <-                                (E/L) - gamma * I
    dR <-                                        gamma * I
    dA <- a * b * c * exp(-(b * exp(-c * times) + c * times))
    return(list(c(dS, dE, dI, dR, dA)))
  })
}

#' SEIRMD Model
#'
#' This function defines the system of ODEs for a SEIRMD model. For use with deSolve::ode.
#'
#' @param times The times at which to evaluate solutions.
#' @param state The current estimate of the variables in the ODE system.
#' @param parameters A vector of parameters. Must include rate parameters \code{beta}, \code{alpha}, \code{gamma}, \code{delta}, \code{lambda}, \code{zeta} and \code{theta}.
#' @export
seirmd <- function(times, state, parameters)
{
  with(as.list(c(state, parameters)), {
    dS <- - (beta * I + alpha * M) * S
    dE <-   (beta * I + alpha * M) * S - lambda * E
    dI <-                                lambda * E - gamma * I - delta * I
    dR <-                                             gamma * I             + zeta * M
    dM <- - theta * M                                           + delta * I - zeta * M
    dD <-   theta * M
    return(list(c(dS, dE, dI, dR, dM, dD)))
  })
}

seirmad <- function(times, state, parameters)
{
    with(as.list(c(state, parameters)), {
        dS <- - (beta1 * I + beta2 * A + alpha * M) * S
        dE <-   (beta1 * I + beta2 * A + alpha * M) * S - lambda1 * E                                                   - lambda2 * E
        dI <-                                            lambda1 * E - gamma1 * I               - delta * I
        dR <-                                                          gamma1 * I  + gamma2 * A             + zeta * M
        dM <- - theta * M                                                                       + delta * I - zeta * M
        dD <-   theta * M
        dA <-                                                                      - gamma2 * A                        + lambda2 * E
        return(list(c(dS, dE, dI, dR, dM, dA, dD)))
    })
}


#' SEIR Model
#'
#' This function defines the system of ODEs for a SEIR model. For use with deSolve::ode.
#'
#' @param times The times at which to evaluate the solutions.
#' @param state The current estimate of the variables in the ODE system.
#' @param parameters A vector of parameters. Must include \code{beta}, \code{rho}, \code{L} and \code{gamma}.
#' @export
seir <- function(times, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- - beta * S * I
    dE <-   beta * S * I - (E/L)
    dI <-                  (E/L) - gamma * I
    dR <-                          gamma * I
    return(list(c(dS, dE, dI, dR)))
  })
}


sir7 <- function(times, state, parameters) # h = 20 for FL, CA; h = 4 for NV
{
  with(as.list(c(state, parameters)), {
    du <- exp( (1/b_s) * exp(b_s * (times + t_s) + a_s) + c_s) * exp(b_s * (times + t_s) + a_s)
    dS <- - du
    dI <-   du - gamma * I
    dR <-   gamma * I
        return(list(c(dS, dI, dR)))
  })
}


