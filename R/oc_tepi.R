#' Compute operating characteristics using TEPI
#'
#' `oc_tepi()` uses the TEPI design to compute operating charateristics of a user-specificed trial scenario.
#' This design maps toxicity and efficacy intervals onto a decision table, forming 16 regions.
#' @param ndose Integer. Number of dose levels. (**Required**)
#' @param target_t Numeric. Target toxicity probability. (**Required**)
#' @param lower_e Numeric. Minimum acceptable efficacy probability. (**Required**)
#' @param ncohort Integer. Number of cohorts. (Default is `10`)
#' @param cohortsize Integer. Size of a cohort. (Default is `3`)
#' @param startdose Integer. Starting dose level. (Default is `1`)
#' @param effint_l Lower efficacy bounds for dose assignment decision table. (Default is `c(0,lower_e,lower_e+0.2,lower_e+0.4)`)
#' @param effint_u Lower efficacy bounds for dose assignment decision table. (Default is `c(lower_e,lower_e+0.2,lower_e+0.4,1)`)
#' @param toxint_l Lower toxicity bounds for dose assignment decision table. (Default is `c(0,0.15,target_t,target_t+0.05)`)
#' @param toxint_u Lower toxicity bounds for dose assignment decision table. (Default is `c(0.15,target_t,target_t+0.05,1)`)
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
#' @examples
#' oc_tepi(
#'   ndose = 5,
#'   target_t = 0.3,
#'   lower_e = 0.4,
#'   ntrial = 10,
#' )
#' @export

# set.seed(30)
oc_tepi <- function(ndose, target_t, lower_e, ncohort = 10, cohortsize = 3,
                    startdose = 1,
                    effint_l = c(0, lower_e, lower_e + 0.2, lower_e + 0.4),
                    effint_u = c(lower_e, lower_e + 0.2, lower_e + 0.4, 1),
                    toxint_l = c(0, 0.15, target_t, target_t + 0.05),
                    toxint_u = c(0.15, target_t, target_t + 0.05, 1),
                    psafe = 0.95, pfutility = 0.95,
                    ntrial = 10000, utilitytype = 1, u1, u2, prob = NULL) {
  safe <- 0
  stop <- 150
  if (utilitytype == 1) {
    u1 <- 60
    u2 <- 40
  }
  if (utilitytype == 2) {
    u1 <- 100
    u2 <- 0
  }



  npts <- ncohort * cohortsize
  YT <- matrix(rep(0, ndose * ntrial), ncol = ndose) # store toxicity outcome
  YE <- matrix(rep(0, ndose * ntrial), ncol = ndose) # store toxicity outcome
  N <- matrix(rep(0, ndose * ntrial), ncol = ndose) # store the number of patients
  dselect <- rep(0, ntrial) # store the selected dose level
  bd.sel <- 0
  bd.pts <- 0
  od.sel <- 0
  od.pts <- 0
  ov.sel <- 0
  ntox <- 0
  neff <- 0
  temp <- get.boundary.utb(target_t, ncohort, cohortsize, cutoff.eli = psafe)
  b.e <- temp[4, ] # escalation boundary
  b.d <- temp[3, ] # deescalation boundary
  b.elim <- temp[2, ] # elimination boundary
  poorall <- 0
  incoherent <- 0
  overdose <- 0
  u.mean <- 0
  tox.l <- seq(0, 0.9, by = 0.1)
  tox.u <- seq(0.1, 1, by = 0.1)
  u.l <- seq(0, 0.9, by = 0.1)
  u.u <- seq(0.1, 1, by = 0.1)
  # effint_l<-c(0,lower_e,lower_e+0.2,lower_e+0.4)
  # effint_u<-c(lower_e,lower_e+0.2,lower_e+0.4,1)
  # toxint_l<-c(0,0.15,target_t,target_t+0.05)
  # toxint_u<-c(0.15,target_t,target_t+0.05,1)
  JUPM <- rep(0, 16)

  ################## simulate trials ###################
  for (trial in 1:ntrial)
  {
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
    pos <- function(pE, yE, yT, n, u1, u2, u) {
      f <- dbeta(pE, yE + 1, n - yE + 1)
      f <- f * pbeta((u1 * pE + u2 - u) / u2, yT + 1, n - yT + 1)
      return(f)
    }

    posH <- rep((u1 * 0.5 + u2) / 10, ndose)
    safe <- 0
    for (i in 1:ncohort)
    {
      ### generate toxicity outcome
      wT <- sum(runif(cohortsize) < pT.true[d])
      yT[d] <- yT[d] + wT
      wE <- sum(runif(cohortsize) < pE.true[d])
      yE[d] <- yE[d] + wE
      n[d] <- n[d] + cohortsize
      nc <- n[d] / cohortsize

      if (!is.na(b.elim[nc])) {
        if (yT[d] >= b.elim[nc]) {
          elimi[d:ndose] <- 1
          if (d == 1) {
            earlystop <- 1
            break
          } else {
            if (sum(elimi[1:d] == 0) > 0) {
              d_opt <- max(which(elimi[1:d] == 0))
              d <- d_opt
              next
            } else {
              earlystop <- 1
              break
            }
          }
        }
      }
      for (ii in 1:4) {
        for (jj in 1:4) {
          a_T <- yT[d] + 1
          b_T <- n[d] - yT[d] + 1
          a_E <- yE[d] + 1
          b_E <- n[d] - yE[d] + 1
          JUPM[(ii - 1) * 4 + jj] <- (pbeta(toxint_u[ii], a_T, b_T) - pbeta(toxint_l[ii], a_T, b_T)) * (pbeta(effint_u[jj], a_E, b_E) - pbeta(effint_l[jj], a_E, b_E))
          JUPM[(ii - 1) * 4 + jj] <- JUPM[(ii - 1) * 4 + jj] / (toxint_u[ii] - toxint_l[ii]) / (effint_u[jj] - effint_l[jj])
        }
      }
      max_JUPM <- which.max(JUPM)
      decision <- 0
      d_opt <- d

      if (max_JUPM <= 7 & d != ndose) {
        if (sum(elimi[(d + 1):ndose] == 0) > 0) {
          decision <- 1
          d_opt <- (d) + min(which(elimi[(d + 1):ndose] == 0))
        }
      }

      if (d > 1) {
        if (max_JUPM == 9 | max_JUPM >= 13) {
          if (sum(elimi[1:(d - 1)] == 0) > 0) {
            decision <- 2
            d_opt <- max(which(elimi[1:d - 1] == 0))
          }
        }
      }

      if (n[d] >= 3 && pbeta(lower_e, yE[d] + 1, n[d] - yE[d] + 1) > pfutility) {
        elimi[d] <- 1
        if (decision == 0) {
          if (sum(elimi[1:d] == 0) > 0) {
            d_opt <- max(which(elimi[1:d] == 0))
          } else {
            earlystop <- 1
            break
          }
        }
        if (decision == 2) {
          if (sum(elimi[1:d] == 0) > 0) {
            d_opt <- max(which(elimi[1:d] == 0))
          } else {
            earlystop <- 1
            break
          }
        }
        if (decision == 1) {
          if (sum(elimi[d:ndose] == 0) > 0) {
            d_opt <- (d - 1) + min(which(elimi[d:ndose] == 0))
          } else if (sum(elimi[1:d] == 0) > 0) {
            d_opt <- max(which(elimi[1:d] == 0))
          } else {
            earlystop <- 1
            break
          }
        }
      }

      if (n[d] >= stop) {
        break
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

      if (abs(u.true[d_opt] - u.true[bd]) <= (0.05 * u.true[bd]) & d_opt <= mtd) {
        od.sel <- od.sel + 1 / ntrial * 100
      }

      if (pT.true[d_opt] > (target_t + 0.1)) {
        ov.sel <- ov.sel + 1 / ntrial * 100
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
