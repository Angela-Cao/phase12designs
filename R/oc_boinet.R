#' Compute Operating Characteristics using BOINET
#'
#' `oc_boinet()` uses the BOINET design to compute operating charateristics of a user-specificed trial scenario.
#' This design uses target toxicity and efficacy rates jointly to form the cutoff intervals within a decision map.
#' @param ndose Integer. Number of dose levels. (**Required**)
#' @param target_t Numeric. Target toxicity probability. (**Required**)
#' @param lower_e Numeric. Minimum acceptable efficacy probability. (**Required**)
#' @param ncohort Integer. Number of cohorts. (Default is `10`)
#' @param cohortsize Integer. Size of a cohort. (Default is `3`)
#' @param startdose Integer. Starting dose level. (Default is `1`)
#' @param psafe Numeric. Early stopping cutoff for toxicity. (Default is `0.95`)
#' @param pfutility Numeric. Early stopping cutoff for efficacy. (Default is `0.95`)
#' @param ntrial  Integer. Number of random trial replications. (Default is `10000`)
#' @param utilitytype Integer. Type of utility structure. (Default is `1`)
#'   - If set to `1`: Use preset weights (w11 = 0.6, w00 = 0.4)
#'   - If set to `2`: Use (w11 = 1, w00 = 0)
#' @param prob Fixed probability vectors. If not specified, a random scenario is used by default.
#' Use this parameter to provide fixed probability vectors as a list with the following named elements:
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
# set.seed(30)
oc_boinet <- function(ndose, target_t, lower_e, ncohort = 10,
                      cohortsize = 3, startdose = 1, psafe = 0.95,
                      pfutility = 0.95, ntrial = 10000, utilitytype = 1,
                      prob = NULL) {
  OBD <- 0
  stop <- 150
  safe <- 0

  if (utilitytype == 1) {
    u1 <- 60
    u2 <- 40
  }
  if (utilitytype == 2) {
    u1 <- 100
    u2 <- 0
  }

  npts <- ncohort * cohortsize
  YT <- matrix(rep(0, ndose * ntrial), ncol = ndose)
  N <- matrix(rep(0, ndose * ntrial), ncol = ndose)
  dselect <- rep(0, ntrial)
  bd.sel <- 0
  bd.pts <- 0
  od.sel <- 0
  od.pts <- 0
  ov.sel <- 0
  ntox <- 0
  neff <- 0
  temp <- get.boundary.utb(target_t, ncohort, cohortsize, cutoff.eli = psafe)
  b.e <- temp[4, ]
  b.d <- temp[3, ]
  b.elim <- temp[2, ]
  poorall <- 0
  incoherent <- 0
  overdose <- 0
  u.mean <- 0

  lambda1 <- 0.16
  lambda2 <- 0.35

  if (target_t == 0.4) {
    lambda1 <- 0.21
    lambda2 <- 0.48
  }
  if (target_t == 0.2) {
    lambda1 <- 0.1
    lambda2 <- 0.23
  }
  eta <- 0.38
  ################## simulate trials ###################
  set.seed(30)

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
    yT <- yE <- rep(0, ndose) ## number of DLT at each dose level
    n <- rep(0, ndose) ## number of patients treated at each dose level
    earlystop <- 0 ## indiate whether the trial terminates early
    d <- startdose ## starting dose level
    elimi <- rep(0, ndose) ## indicate whether doses are eliminated due to toxicity
    elimiE <- rep(0, ndose) ## indicate whether doses are eliminated due to efficacy
    incoh <- 0 ## count incoherent movement


    posH <- rep((u1 * 0.5 + u2) / 10, ndose)
    safe <- 0
    for (i in 1:ncohort) {
      wT <- sum(runif(cohortsize) < pT.true[d])
      yT[d] <- yT[d] + wT
      wE <- sum(runif(cohortsize) < pE.true[d])
      yE[d] <- yE[d] + wE
      n[d] <- n[d] + cohortsize
      nc <- n[d] / cohortsize
      if (n[d] >= stop) {
        break
      }
      if (!is.na(b.elim[nc])) {
        if (yT[d] >= b.elim[nc]) {
          elimi[d:ndose] <- 1
          if (d == 1) {
            earlystop <- 1
            break
          }
        }
      }
      if (n[d] >= 3 && pbeta(lower_e, yE[d] + 1, n[d] - yE[d] + 1) > pfutility) {
        elimi[d] <- 1
      }
      phatT <- yT / (n + 0.0000001) + runif(ndose) * 10^(-10)
      phatE <- yE / (n + 0.0000001) + runif(ndose) * 10^(-10)
      phatE <- phatE * (1 - elimi)
      if (phatT[d] >= lambda2 && d != 1) {
        if (sum(elimi[1:(d - 1)] == 0) > 0) {
          d_opt <- max(which(elimi[1:(d - 1)] == 0))
        } else {
          d_opt <- d
        }
      } else if (phatT[d] >= lambda2 && d == 1) {
        if (elimi[d] == 0) {
          d_opt <- d
        } else {
          earlystop <- 1
          break
        }
      } else {
        if (phatE[d] > eta) {
          d_opt <- d
        } else if (phatT[d] <= lambda1) {
          if (sum(elimi[min(d + 1, ndose):ndose] == 0) > 0) {
            d_opt <- min(d + min(which(elimi[min(d + 1, ndose):ndose] == 0)), ndose)
          } else {
            d_opt <- d
          }
        } else {
          admi_set <- d
          if (d > 1) {
            if (sum(elimi[1:(d - 1)] == 0) > 0) {
              admi_set <- c(admi_set, max(which(elimi[1:(d - 1)] == 0)))
            }
          }
          if (d < ndose) {
            if (sum(elimi[(d + 1):ndose] == 0) > 0) {
              admi_set <- c(admi_set, d + min(which(elimi[(d + 1):ndose] == 0)))
            }
          }
          if (n[admi_set[length(admi_set)]] == 0) {
            d_opt <- admi_set[length(admi_set)]
          } else {
            temp_phat <- phatE[admi_set]
            d_opt <- admi_set[which.max(temp_phat)]
          }
        }
      }


      if (elimi[d_opt] == 1) {
        earlystop <- 1
        break
      }
      if (sum(elimi) == ndose) {
        earlystop <- 1
        break
      }
      if (((yT[d] / n[d]) > target_t) & d_opt > d) {
        incoh <- incoh + 1
      }
      d <- d_opt
    }
    incoherent <- incoherent + (incoh / i) / ntrial * 100
    if (earlystop == 0) {
      pT <- (yT + 0.05) / (n + 0.1)
      pE <- (yE + 0.05) / (n + 0.1)
      pT <- pava(pT, n + 0.1) + 0.001 * seq(1, ndose)
      pE <- peestimate(yE, n)
      u <- u1 * pE + (1 - pT) * u2
      u[elimi == 1] <- -100
      u[elimiE == 1] <- -100
      u[n == 0] <- -100
      # u[pT>(target_t+0.1)]<--100
      d_mtd <- which.min(abs(pT - target_t))
      d_opt <- which.max(u[1:d_mtd])
      dselect[trial] <- d_opt
      if (d_opt == bd) {
        bd.sel <- bd.sel + 1 / ntrial * 100
      }
      if (pT.true[d_opt] > (target_t + 0.1)) {
        ov.sel <- ov.sel + 1 / ntrial * 100
      }
      if (abs(u.true[d_opt] - u.true[bd]) <= (0.05 * u.true[bd]) & d_opt <= mtd) {
        od.sel <- od.sel + 1 / ntrial * 100
      }
      # u.mean<-u.mean+utility[d_opt]/ntrial
    } else {
      dselect[trial] <- 99
    }
    earlystop <- sum(dselect == 99) / ntrial * 100
    if (n[bd] < (npts / ndose)) {
      poorall <- poorall + 1 / ntrial * 100
    }
    overdose <- overdose + sum(n[pT.true > (target_t + 0.1)]) / ntrial / npts * 100
    bd.pts <- bd.pts + n[bd] / ntrial / npts * 100
    od.pts <- od.pts + sum(n[abs(u.true[1:mtd] - (u.true[bd])) <= (0.05 * u.true[bd])]) / ntrial / npts * 100
    ntox <- ntox + sum(yT) / ntrial
    neff <- neff + sum(yE) / ntrial
  }
  results <- list(
    bd.sel = bd.sel, od.sel = od.sel, bd.pts = bd.pts, od.pts = od.pts,
    earlystop = earlystop, ntox = ntox, neff = neff,
    overdose = overdose, poorall = poorall, incoherent = incoherent, ov.sel = ov.sel
  )
  return(results)
}
