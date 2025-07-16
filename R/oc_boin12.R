#' Compute Operating Characteristics using BOIN12
#'
#' `oc_boin12()` uses the BOIN12 design to compute operating charateristics of a user-specificed trial scenario.
#' This design places significance on optimizing utility and the toxicity–efficacy trade-off.
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
#' ``
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
oc_boin12 <- function(ndose, target_t, lower_e, ncohort = 10,
                      cohortsize = 3, startdose = 1,
                      psafe = 0.95, pfutility = 0.95,
                      ntrial = 10000, utilitytype = 1,
                      u1, u2, prob = NULL) {
  if (utilitytype == 1) {
    u1 <- 60
    u2 <- 40
  } else if (utilitytype == 2) {
    u1 <- 100
    u2 <- 0
  }
  n.earlystop <- ncohort * cohortsize
  p.saf <- 0.6 * target_t
  p.tox <- 1.4 * target_t
  N1 <- 6
  N2 <- 9
  randomtype <- 1
  poorall <- 0
  incoherent <- 0
  overdose <- 0
  bd.sel <- 0
  bd.pts <- 0
  od.sel <- 0
  ov.sel <- 0
  od.pts <- 0
  ntox <- 0
  neff <- 0
  npts <- ncohort * cohortsize
  YT <- matrix(rep(0, ndose * ntrial), ncol = ndose) # store toxicity outcome
  YE <- matrix(rep(0, ndose * ntrial), ncol = ndose) # store efficacy outcome
  N <- matrix(rep(0, ndose * ntrial), ncol = ndose) # store the number of patients
  dselect <- rep(0, ntrial) # store the selected dose level
  durationV <- rep(0, ntrial)
  sel <- rep(0, ndose)
  pts <- rep(0, ndose)
  dlt <- rep(0, ndose)
  eff <- rep(0, ndose)
  poorall <- 0
  temp <- get.boundary(target_t, lower_e, ncohort, cohortsize, cutoff.eli = psafe, cutoff.eli.E = pfutility)
  b.e <- temp[4, ] # escalation boundary
  b.d <- temp[3, ] # deescalation boundary
  b.elim <- temp[2, ] # elimination boundary
  b.elimE <- temp[5, ]
  u01 <- 100
  u10 <- 0
  u11 <- u1
  u00 <- u2
  utility <- c(u11, u10, u01, u00)
  # Assume independence between toxicity and efficacy
  targetP <- c(lower_e * target_t, target_t * (1 - lower_e), (1 - target_t) * lower_e, (1 - target_t) * (1 - lower_e))
  # Calculate the benchmark utility
  uu <- sum(targetP * utility) # highest unacceptable utility
  uu <- uu + (100 - uu) / 2 # benchmark utility (i.e., desirable utility)
  # Calculate true utility
  p10 <- p01 <- p00 <- p11 <- rep(0, ndose)
  if (FALSE) {
    for (d in 1:ndose) {
      p11[d] <- integrate(f1, bn.m1 = qnorm(pT.true[d]), bn.m2 = qnorm(pE.true[d]), rho = rho, lower = 0, upper = Inf)$value
      p10[d] <- pT.true[d] - p11[d]
      p01[d] <- pE.true[d] - p11[d]
      p00[d] <- 1 - p11[d] - p10[d] - p01[d]
    }
  }
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
    yT <- yE <- rep(0, ndose) ## number of DLT/efficacy at each dose level
    y01 <- y10 <- y11 <- y00 <- rep(0, ndose) ## number of different outcomes at each dose level
    n <- rep(0, ndose) ## number of patients treated at each dose level
    earlystop <- 0 ## indicate whether the trial terminates early
    d <- startdose ## starting dose level
    elimi <- rep(0, ndose) ## whether doses are eliminated due to toxicity
    elimiE <- rep(0, ndose) ## whether doses are eliminated due to efficacy
    safe <- 0
    posH <- rep(1 - uu / 100, ndose)
    duration <- 0
    for (i in 1:ncohort) {
      T.time <- 0 # obscohort$t.tox
      E.time <- 0 # obscohort$t.tox
      inter.arrival <- 0 # cumsum(rexp(cohortsize,rate=accrual.rate))
      t.all.seen <- 0 # inter.arrival+pmax(T.time,E.time)
      duration <- duration + max(t.all.seen)
      n[d] <- n[d] + cohortsize
      wT <- sum(runif(cohortsize) < pT.true[d])
      yT[d] <- yT[d] + wT
      wE <- sum(runif(cohortsize) < pE.true[d])
      yE[d] <- yE[d] + wE
      nc <- n[d] / cohortsize
      # determine whether current dose level is overly toxic
      if (!is.na(b.elim[nc])) {
        if (yT[d] >= b.elim[nc]) {
          elimi[d:ndose] <- 1
          if (d == 1) {
            earlystop <- 1
            break
          }
        }
      }
      if (!is.na(b.elimE[nc])) {
        if (yE[d] <= b.elimE[nc]) {
          elimi[d] <- 1
        }
      }
      if (sum(elimi == 1) == ndose) {
        earlystop <- 1
        break
      }
      u_curr <- (u1 * yE[d] / n[d] + (1 - yT[d] / n[d]) * u2) / 100 * n[d]
      posH[d] <- 1 - pbeta(uu / 100, 1 + u_curr, n[d] - u_curr + 1)
      posH <- posH * (1 - elimi)
      if (n[d] >= N1) {
        safe <- 1
      } else {
        safe <- 0
      }
      if (n[d] >= n.earlystop) {
        break
      }
      if (yT[d] >= b.d[nc] && d != 1) {
        if (sum(elimi[1:(d - 1)] == 0) > 0) {
          d_opt <- max(which(elimi[1:(d - 1)] == 0))
        } else {
          if (elimi[d] == 1) {
            earlystop <- 1
            break
          } else {
            d_opt <- d
          }
        }
      } else if (yT[d] >= b.d[nc] && d == 1) {
        if (elimi[d] == 0) {
          d_opt <- d
        } else {
          earlystop <- 1
          break
        }
      } else {
        admi_set <- d
        if (d > 1) {
          if (sum(elimi[1:(d - 1)] == 0) > 0) {
            admi_set <- c(admi_set, max(which(elimi[1:(d - 1)] == 0)))
          }
        }
        if (d < ndose) {
          if (safe == 0) {
            if (sum(elimi[(d + 1):ndose] == 0) > 0) {
              admi_set <- c(admi_set, d + min(which(elimi[(d + 1):ndose] == 0)))
            }
          } else {
            if (yT[d] <= b.e[nc] && sum(elimi[(d + 1):ndose] == 0) > 0) {
              admi_set <- c(admi_set, d + min(which(elimi[(d + 1):ndose] == 0)))
            }
          }
        }
        temp.posH <- posH[admi_set] + runif(length(admi_set)) * (10^-15)
        d_opt <- admi_set[which.max(temp.posH)]
      }

      if (elimi[d_opt] == 1) {
        earlystop <- 1
        break
      }
      if (sum(elimi) == ndose) {
        earlystop <- 1
        break
      }
      if (d < ndose) {
        if (sum(elimi[(d + 1):ndose] == 0) > 0) {
          d_temp <- d + min(which(elimi[(d + 1):ndose] == 0))
          if (n[d] >= N2 && n[min(d_temp, ndose)] == 0 && yT[d] < b.d[n[d] / cohortsize]) { # nolint: line_length_linter.
            d_opt <- d_temp
          }
        }
      }
      d <- d_opt
    }
    YT[trial, ] <- yT # nolint: object_name_linter.
    YE[trial, ] <- yE
    N[trial, ] <- n
    durationV[trial] <- duration
    if (earlystop == 0) {
      pT <- (yT + 0.05) / (n + 0.1)
      pE <- (yE + 0.05) / (n + 0.1)
      pT <- pava(pT, n + 0.1) + 0.001 * seq(1, ndose)
      pE <- peestimate(yE, n)
      u <- u1 * pE + (1 - pT) * u2
      u[elimi == 1] <- -100
      u[elimiE == 1] <- -100
      u[n == 0] <- -100
      d_mtd <- which.min(abs(pT - target_t))
      d_opt <- which.max(u[1:d_mtd])
      dselect[trial] <- d_opt
      if (d_opt == bd) {
        bd.sel <- bd.sel + 1 / ntrial * 100
      }
      if (abs(u.true[d_opt] - u.true[bd]) <= (0.05 * u.true[bd]) && d_opt <= mtd) {
        od.sel <- od.sel + 1 / ntrial * 100
      }
      if (pT.true[d_opt] > (target_t + 0.1)) {
        ov.sel <- ov.sel + 1 / ntrial * 100
      }
      sel[dselect[trial]] <- sel[dselect[trial]] + 1 / ntrial * 100
    } else {
      dselect[trial] <- 99
    }
    earlystop <- sum(dselect == 99) / ntrial * 100
    if (n[bd] < (npts / ndose)) {
      poorall <- poorall + 1 / ntrial * 100
    }
    overdose <- overdose + sum(n[pT.true > (target_t + 0.1)]) / ntrial / npts * 100
    bd.pts <- bd.pts + n[bd] / ntrial / npts * 100
    od.pts <- od.pts + sum(n[abs(u.true[1:mtd] - u.true[bd]) <= (0.05 * u.true[bd])]) / ntrial / npts * 100
    pts <- pts + n / ntrial
    dlt <- dlt + yT / ntrial
    eff <- eff + yE / ntrial
    ntox <- ntox + sum(yT) / ntrial
    neff <- neff + sum(yE) / ntrial
  }
  sel <- round(sel, 1)
  pts <- round(pts, 1)
  dlt <- round(dlt, 1)
  u.true <- round(u.true, 1)
  earlystop <- sum(dselect == 99) / ntrial * 100
  results <- list(
    bd.sel = bd.sel, od.sel = od.sel, bd.pts = bd.pts, od.pts = od.pts,
    earlystop = earlystop, ntox = ntox, neff = neff, u.mean = 0,
    overdose = overdose, poorall = poorall, incoherent = 0, ov.sel = ov.sel
  )
  results
}
