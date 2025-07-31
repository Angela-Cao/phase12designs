#' Compute operating characteristics using STEIN
#'
#' `oc_stein()` uses the STEIN design to compute operating charateristics of a user-specificed trial scenario.
#' This design uses target toxicity and efficacy rates separately to form the cutoff intervals within a decision map.
#' @param ndose Integer. Number of dose levels. (**Required**)
#' @param target_t Numeric. Target toxicity probability. (**Required**)
#' @param lower_e Numeric. Minimum acceptable efficacy probability. (**Required**)
#' @param ncohort Integer. Number of cohorts. (Default is `10`)
#' @param cohortsize Integer. Size of a cohort. (Default is `3`)
#' @param startdose Integer. Starting dose level. (Default is `1`)
#' @param psi1 Numerical. Highest inefficacious efficacy probability.
#' @param psi2 Numerical. Lowest highly-promising efficacy probability.
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
#' oc_stein(
#'   ndose = 5,
#'   target_t = 0.3,
#'   lower_e = 0.4,
#'   ntrial = 10,
#' )
#' @export
oc_stein <- function(ndose, target_t, lower_e,
                     ncohort = 10, cohortsize = 3, startdose = 1,
                     psi1 = 0.2, psi2 = 0.6, psafe = 0.95, pfutility = 0.9,
                     ntrial = 10000, utilitytype = 1, u1, u2,
                     prob = NULL) {
  if (utilitytype == 1) {
    u1 <- 60
    u2 <- 40
  }
  if (utilitytype == 2) {
    u1 <- 100
    u2 <- 0
  }

  npts <- ncohort * cohortsize
  YT <- matrix(0, ncol = ndose, nrow = ntrial)
  YE <- matrix(0, ncol = ndose, nrow = ntrial)
  N <- matrix(0, ncol = ndose, nrow = ntrial)
  dselect <- rep(0, ntrial)

  sel <- pts <- dlt <- eff <- rep(0, ndose)

  ntox <- neff <- acr <- exc <- 0
  temp <- stein.boundary(target_t, ncohort, cohortsize)
  b.e <- temp[4, ]
  b.d <- temp[3, ]
  b.elim <- temp[2, ]

  psi <- log((1 - psi1) / (1 - psi2)) / log(psi2 * (1 - psi1) / (psi1 * (1 - psi2)))

  bd.sel <- bd.pts <- od.sel <- od.pts <- ov.sel <- 0
  poorall <- incoherent <- overdose <- 0
  # set.seed(30)
  # Simulate trials
  for (trial in 1:ntrial) {
    yT <- yE <- n <- rep(0, ndose)
    earlystop <- 0
    d <- startdose
    elimi <- elimiE <- rep(0, ndose)

    if (!is.null(prob)) {
      probs <- prob
    } else {
      probs <- simprob(ndose, lower_e, target_t, u1, u2, randomtype)
    }
    pE.true <- probs$pE
    pT.true <- probs$pT
    u.true <- u1 * pE.true + (1 - pT.true) * u2
    bd <- probs$obd
    mtd <- probs$mtd

    for (i in seq_len(ncohort)) {
      wT <- sum(runif(cohortsize) < pT.true[d])
      wE <- sum(runif(cohortsize) < pE.true[d])

      yT[d] <- yT[d] + wT
      yE[d] <- yE[d] + wE
      n[d] <- n[d] + cohortsize
      nc <- n[d] / cohortsize

      if (!is.na(b.elim[nc]) && yT[d] >= b.elim[nc]) {
        elimi[d:ndose] <- 1
        if (d == 1) {
          earlystop <- 1
          break
        }
      }

      if (n[d] >= 3 && pbeta(lower_e, yE[d] + 1, n[d] - yE[d] + 1) > pfutility) {
        elimiE[d] <- 1
      }

      if (yT[d] >= b.d[nc] && d != 1) {
        d_opt <- max(which(elimi[1:(d - 1)] == 0))
      } else if (yT[d] >= b.d[nc] && d == 1) {
        d_opt <- d
      } else {
        if (yE[d] / n[d] > psi) {
          d_opt <- d
        } else if (yT[d] > b.e[nc]) {
          posH <- c(
            pbeta(psi, yE[max(1, d - 1)] + 1, n[max(1, d - 1)] - yE[max(1, d - 1)] + 1),
            pbeta(psi, yE[d] + 1, n[d] - yE[d] + 1)
          ) - c(0, 0.001)
          d_opt <- max(d - 2 + which.min(posH), 1)
        } else if (yT[d] <= b.e[nc] && d == 1) {
          if (elimi[d + 1] == 0) {
            posH <- c(
              pbeta(psi, yE[d] + 1, n[d] - yE[d] + 1),
              ifelse(n[d + 1] > 0, pbeta(psi, yE[d + 1] + 1, n[d + 1] - yE[d + 1] + 1), 0)
            ) - c(0, 0.001)
            d_opt <- d - 1 + which.min(posH)
          } else {
            d_opt <- d
          }
        } else if (yT[d] <= b.e[nc] && d == ndose) {
          posH <- c(
            pbeta(psi, yE[d - 1] + 1, n[d - 1] - yE[d - 1] + 1),
            pbeta(psi, yE[d] + 1, n[d] - yE[d] + 1)
          ) - c(0, 0.001)
          d_opt <- d - 2 + which.min(posH)
        } else if (yT[d] <= b.e[nc]) {
          if (elimi[d + 1] == 0) {
            posH <- c(
              pbeta(psi, yE[d - 1] + 1, n[d - 1] - yE[d - 1] + 1),
              pbeta(psi, yE[d] + 1, n[d] - yE[d] + 1),
              ifelse(n[d + 1] > 0, pbeta(psi, yE[d + 1] + 1, n[d + 1] - yE[d + 1] + 1), 0)
            ) - c(0, 0.001, 0.002)
            d_opt <- d - 2 + which.min(posH)
          } else {
            posH <- c(
              pbeta(psi, yE[d - 1] + 1, n[d - 1] - yE[d - 1] + 1),
              pbeta(psi, yE[d] + 1, n[d] - yE[d] + 1)
            ) - c(0, 0.001)
            d_opt <- d - 2 + which.min(posH)
          }
        }
      }

      if (elimiE[d_opt] == 1) {
        earlystop <- 1
        break
      }

      d <- d_opt
    }

    YT[trial, ] <- yT
    YE[trial, ] <- yE
    N[trial, ] <- n

    ntox <- ntox + sum(yT) / ntrial
    neff <- neff + sum(yE) / ntrial

    if (earlystop == 0) {
      pT <- (yT + 0.05) / (n + 0.1)
      pE <- (yE + 0.05) / (n + 0.1)
      pT <- pava(pT, n + 0.1) + 0.001 * seq_len(ndose)
      pE <- peestimate(yE, n)

      u <- u1 * pE + (1 - pT) * u2
      u[elimi == 1 | elimiE == 1 | n == 0] <- -100

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

    pts <- pts + n / ntrial
    dlt <- dlt + yT / ntrial
    eff <- eff + yE / ntrial

    if (n[bd] < (npts / ndose)) {
      poorall <- poorall + 1 / ntrial * 100
    }
    overdose <- overdose + sum(n[pT.true > (target_t + 0.1)]) / ntrial / npts * 100
    bd.pts <- bd.pts + n[bd] / ntrial / npts * 100
    od.pts <- od.pts + sum(n[abs(u.true[1:mtd] - u.true[bd]) <= (0.05 * u.true[bd])]) / ntrial / npts * 100
  }

  pts <- round(pts, 1)
  dlt <- round(dlt, 1)

  results <- list(
    bd.sel = bd.sel, od.sel = od.sel,
    bd.pts = bd.pts, od.pts = od.pts,
    earlystop = sum(dselect == 99) / ntrial * 100,
    ntox = ntox, neff = neff, u.mean = 0,
    overdose = overdose, poorall = poorall,
    incoherent = 0, ov.sel = ov.sel
  )

  return(results)
}
