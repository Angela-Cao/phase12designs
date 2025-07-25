#' Compute operating characteristics using Ji3+3
#'
#' `oc_ji3p3()` uses the  Ji3+3 design to compute operating charateristics of a user-specificed trial scenario.
#' This design compares observed efficacy and toxicity with predefined target rates.
#' @param ndose Integer. Number of dose levels. (**Required**)
#' @param target_t Numeric. Target toxicity probability. (**Required**)
#' @param target_e Numeric. Target efficacy probability. (**Required**)
#' @param lower_e Numeric. Minimum acceptable efficacy probability. (**Required**)
#' @param ncohort Integer. Number of cohorts. (Default is `10`)
#' @param cohortsize Integer. Size of a cohort. (Default is `3`)
#' @param startdose Integer. Starting dose level. (Default is `1`)
#' @param eps1 Numerical. Width of the subrectangle.
#' @param eps2 Numerical. Width of the subreactangle.
#' @param psafe Numeric. Early stopping cutoff for toxicity. (Default is `0.95`)
#' @param pfutility Numeric. Early stopping cutoff for efficacy. (Default is `0.95`)
#' @param ntrial  Integer. Number of random trial replications. (Default is `10000`)
#' @param utilitytype Integer. Type of utility structure. (Default is `1`)
#'   - If set to `1`: Use preset weights (w11 = 0.6, w00 = 0.4)
#'   - If set to `2`: Use (w11 = 1, w00 = 0)
#'   - Other: Use user-specified values from `u1` and `u2`.
#' @param u1 Numeric. Utility parameter w_11. (0-100)
#' @param u2 Numeric. Utility parameter w_00. (0-100)
#' @param prob Fixed probability vectors. If not specified, a random scenario is used by default.
#' Use this parameter to provide fixed probability vectors as a list of the following named elements:
#'   - `pE`: Numeric vector of efficacy probabilities for each dose level.
#'   - `pT`: Numeric vector of toxicity probabilities for each dose level.
#'   - `obd`: Integer indicating the index of the true Optimal Biological Dose (OBD).
#'   - `mtd`: Integer indicating the index of the true Maximum Tolerated Dose (MTD).
#'
#' For example:
#' ```r
#' prob <- list(
#'   pE = c(0.4, 0.5, 0.6, 0.6, 0.6),
#'   pT = c(0.1, 0.2, 0.3, 0.4, 0.4),
#'   obd = 3,
#'   mtd = 2
#' )
#' ```
#' @return A list containing operating characteristics such as:
#' \describe{
#'   \item{bd.sel}{OBD selection percentage}
#'   \item{od.sel}{Favorable dose selection percentage}
#'   \item{bd.pts}{Average percentage of patients at the OBD }
#'   \item{od.pts}{Average percentage of patients at the favorable doses}
#'   \item{earlystop}{Percentage of early stopped trials}
#'   \item{overdose}{Overdose patients percentage }
#'   \item{poorall}{Poor allocation percentage}
#'   \item{ov.sel}{Overdose selection percentage}
#' }
#'
#' @export

oc_ji3p3 <- function(ndose, target_t, target_e,
                     lower_e = 0.2, ncohort = 10, cohortsize = 3, startdose = 1,
                     eps1 = 0.05, eps2 = 0.05,
                     psafe = 0.95, pfutility = 0.95,
                     ntrial = 10000, utilitytype = 1, u1, u2,
                     prob = NULL) {
  if (utilitytype == 1) {
    u1 <- 60
    u2 <- 40
  } else if (utilitytype == 2) {
    u1 <- 100
    u2 <- 0
  }
  parabeta <- c(1, 1) # Beta distribution parameters
  selection <- numeric(ntrial)
  pat_treated <- matrix(NA, nrow = ntrial, ncol = ndose)
  csize <- cohortsize
  samplesize <- csize * ncohort
  targetT <- target_t
  targetE <- lower_e
  npts <- ncohort * cohortsize
  bd.sel <- 0
  bd.pts <- 0
  od.sel <- 0
  od.pts <- 0
  ov.sel <- 0
  ntox <- 0
  neff <- 0
  poorall <- 0
  incoherent <- 0
  overdose <- 0
  dselect <- rep(0, ntrial) # store the selected dose level

  for (trial in 1:ntrial) {
    if (!is.null(prob)) {
      probs <- prob
    } else {
      probs <- simprob(ndose, lower_e, target_t, u1, u2, randomtype)
    }
    jj <- probs$pE
    kk <- probs$pT
    pE.true <- jj
    pT.true <- kk
    u.true <- (u1 * pE.true + (1 - pT.true) * u2)
    bd <- probs$obd
    mtd <- probs$mtd
    toxtrue <- pT.true
    efftrue <- pE.true
    x <- rep(0, ndose) # toxicity number for each dose
    y <- rep(0, ndose) # efficacy number for each dose
    n <- rep(0, ndose) # patient number for each dose
    d <- startdose
    st <- 0 # stop sign
    currentsize <- 0
    dtox <- rep(0, ndose)
    dfut <- rep(0, ndose)
    earlystop <- 0
    while (st == 0) {
      xx <- sum(runif(cohortsize) < toxtrue[d])
      yy <- sum(runif(cohortsize) < efftrue[d])
      x[d] <- x[d] + xx # total tox outcome
      y[d] <- y[d] + yy # total eff outcome
      n[d] <- n[d] + csize # total sample size used
      currentsize <- currentsize + csize
      # safety rule	(monotonic tox)
      if (pbeta(target_t, parabeta[1] + x[d], parabeta[2] + n[d] - x[d]) < 1 - psafe) {
        dtox[d:ndose] <- 1
      }
      # futility rule (mark only current dose)
      if (pbeta(lower_e, parabeta[1] + y[d], parabeta[2] + n[d] - y[d]) > pfutility) {
        dfut[d] <- 2
      }
      # obtain optimal dosing decisions
      dec <- ji3plus3(x = x[d], y = y[d], n = n[d], pT = target_t, eps1 = eps1, eps2 = eps2, pE = target_e)
      # if the sampsize is reached, stop
      if (currentsize >= samplesize) {
        st <- 1
      } else {
        # if all doses are either too toxic or of no efficacy (couldn't find a dose for the next corhort)
        if (length(which(dfut + dtox == 0)) == 0) {
          st <- 1
          earlystop <- 1
        } else {
          # if current dose is too toxic
          if (dtox[d] == 1) {
            # if there is a lower dose available, de-escalate to the closet dose below
            if (d != 1 && length(which(dtox[1:(d - 1)] + dfut[1:(d - 1)] == 0) != 0)) {
              d <- max(which(dtox[1:(d - 1)] + dfut[1:(d - 1)] == 0))
            } else {
              st <- 1 # if there is no valid dose below current dose
            }
          }
          # if current dose is of no efficacy
          else if (dfut[d] == 2) {
            if (dec == "E") {
              # if there is a valid dose above current dose, escalate by minimum dose size
              if (d != ndose && length(which(dtox[(d + 1):ndose] + dfut[(d + 1):ndose] == 0)) != 0) {
                d <- min(which(dtox[(d + 1):ndose] + dfut[(d + 1):ndose] == 0)) + d
              } # if there is no dose above current dose, de-escalate to the next available dose
              else if (d != 1 && length(which(dtox[1:(d - 1)] + dfut[1:(d - 1)] == 0) != 0)) {
                d <- max(which(dtox[1:(d - 1)] + dfut[1:(d - 1)] == 0))
              } # if there is no dose above or below current dose, terminate the trial
              else {
                st <- 1
                earlystop <- 1
              }
            } else if (dec == "D") {
              # if there is a valid dose below current dose, de-escalate by minimum dose size
              if (d != 1 && length(which(dtox[1:(d - 1)] + dfut[1:(d - 1)] == 0) != 0)) {
                d <- max(which(dtox[1:(d - 1)] + dfut[1:(d - 1)] == 0))
              } else {
                st <- 1 # if there is no valid dose below current dose
                earlystop <- 1
              }
            } else {
              # when dec = "S", if there are valid doses above the current dose, escalate
              # if not, de-escalate, if neither, stop
              if (sum(which(dfut + dtox == 0) < d) != 0) {
                d <- max(which(dtox[1:(d - 1)] + dfut[1:(d - 1)] == 0))
              } else {
                st <- 1
                earlystop <- 1
              }
            }
          }
          # Below is the condition where the current dose is efficacious and not too toxic
          else {
            if (dec == "E") {
              # if there is a higher dose available, escalate
              if (d != ndose && length(which(dtox[(d + 1):ndose] + dfut[(d + 1):ndose] == 0)) != 0) {
                d <- min(which(dtox[(d + 1):ndose] + dfut[(d + 1):ndose] == 0)) + d
              } else {
                d <- d
              }
            } else if (dec == "D") {
              # if there is a lower dose available, de-escalate
              if (d != 1 && (length(which(dtox[1:(d - 1)] + dfut[1:(d - 1)] == 0) != 0))) {
                d <- max(which(dtox[1:(d - 1)] + dfut[1:(d - 1)] == 0))
              } # else?
              else {
                d <- d
              }
            } else if (dec == "S") {
              d <- d
            }
          }
        }
      }
    }

    if (earlystop == 0) {
      yT_c <- x
      yE_c <- y
      elimi <- dtox
      elimiE <- ifelse(dfut == 2, 1, 0)
      pT_est <- (yT_c + 0.05) / (n + 0.1)
      pE_est <- (yE_c + 0.05) / (n + 0.1)
      pT_est <- pava(pT_est, n + 0.1) + 0.001 * seq(1, ndose)
      pE_est <- peestimate(yE_c, n)
      u <- u1 * pE_est + (1 - pT_est) * u2
      u[elimi == 1] <- -100
      u[elimiE == 1] <- -100
      u[n == 0] <- -100
      d_mtd <- which.min(abs(pT_est - targetT))
      d_opT_est <- which.max(u[1:d_mtd])
      dselect[trial] <- d_opT_est
      if (d_opT_est == bd) {
        bd.sel <- bd.sel + 1 / ntrial * 100
      }
      if (abs(u.true[d_opT_est] - u.true[bd]) <= (0.05 * u.true[bd]) & d_opT_est <= mtd) {
        od.sel <- od.sel + 1 / ntrial * 100
      }
      if (pT.true[d_opT_est] > (targetT + 0.1)) {
        ov.sel <- ov.sel + 1 / ntrial * 100
      }
      dselect[trial] <- d_opT_est
    } else {
      dselect[trial] <- 99
    }

    earlystop <- sum(dselect == 99) / ntrial * 100
    if (n[bd] < (npts / ndose)) {
      poorall <- poorall + 1 / ntrial * 100
    }
    overdose <- overdose + sum(n[pT.true > (targetT + 0.1)]) / ntrial / npts * 100
    bd.pts <- bd.pts + n[bd] / ntrial / npts * 100
    od.pts <- od.pts + sum(n[abs(u.true[1:mtd] - u.true[bd]) <= (0.05 * u.true[bd])]) / ntrial / npts * 100
  }


  results <- list(
    bd.sel = bd.sel, od.sel = od.sel, bd.pts = bd.pts, od.pts = od.pts,
    earlystop = earlystop, ntox = 0, neff = 0, u.mean = 0,
    overdose = overdose, poorall = poorall, incoherent = 0, ov.sel = ov.sel
  )

  return(results)
}
