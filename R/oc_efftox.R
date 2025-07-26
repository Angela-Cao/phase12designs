#' Compute operating characteristics using EffTox
#'
#' `oc_efftox()` uses the EffTox design to compute operating charateristics of a user-specificed trial scenario.
#' This design uses toxicity–efficacy trade-off contours.
#' @param ndose Integer. Number of dose levels. (**Required**)
#' @param target_t Numeric. Target toxicity probability. (**Required**)
#' @param lower_e Numeric. Minimum acceptable efficacy probability. (**Required**)
#' @param ncohort Integer. Number of cohorts. (Default is `10`)
#' @param startdose Integer. Starting dose level. (Default is `1`)
#' @param ntrial  Integer. Number of random trial replications. (Default is `10000`)
#' @param utilitytype Integer. Type of utility structure. (Default is `1`)
#'   - If set to `1`: Use preset weights (w11 = 0.6, w00 = 0.4)
#'   - If set to `2`: Use (w11 = 1, w00 = 0)
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

oc_efftox <- function(ndose, target_t, lower_e,
                      ncohort = 10, startdose = 1, ntrial = 10000,
                      utilitytype = 1, prob = NULL) {
  if (utilitytype == 1) {
    u1 <- 60
    u2 <- 40
    eff0 <- 0.20
    tox1 <- 0.533
    eff_star <- 0.5
    tox_star <- 0.20
  } else if (utilitytype == 2) {
    u1 <- 100
    u2 <- 0
    eff0 <- 0.80
    tox1 <- 0.99
    eff_star <- 0.90
    tox_star <- 0.50
  }

  targetT <- target_t
  targetE <- lower_e
  cohortsize <- 3
  npts <- cohortsize * ncohort
  uu <- u1 * 0.7 + (1 - 0.1) * u2
  cutoff.eliT <- 0.95
  cutoff.eliE <- 0.9
  res <- NULL
  dselect <- rep(0, ntrial) # store the selected dose level
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
  u.mean <- 0
  pp <- efftox_solve_p(
    eff0 = eff0, tox1 = tox1,
    eff_star = eff_star, tox_star = tox_star
  )
  dat <- list(
    num_doses = ndose, real_doses = seq(1, ndose, by = 1),
    efficacy_hurdle = lower_e, toxicity_hurdle = target_t,
    p_e = 0.05, p_t = 0.1, p = pp,
    eff0 = eff0, tox1 = tox1, eff_star = eff_star, tox_star = tox_star,
    alpha_mean = -2.823, alpha_sd = 2.7099,
    beta_mean = 3.9364, beta_sd = 2.7043,
    gamma_mean = -2.8240, gamma_sd = 2.7108,
    zeta_mean = 3.9374, zeta_sd = 2.7032,
    eta_mean = 0, eta_sd = 0.2,
    psi_mean = 0, psi_sd = 1,
    doses = c(),
    tox = c(),
    eff = c(),
    num_patients = 0
  )

  for (trial in 1:ntrial) {
    set.seed(30 + trial)
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
    temp <- efftox_simulate(dat,
      num_sims = 1, first_dose = startdose,
      true_eff = pE.true, true_tox = pT.true,
      cohort_sizes = rep(cohortsize, ncohort),
      chains = 1, iter = 500, show_messages = FALSE,
      open_progress = FALSE, refresh = 0
    )

    if (is.na(temp$recommended_dose)) {
      dselect[trial] <- 99
    } else {
      dselect[trial] <- temp$recommended_dose
    }

    if (dselect[trial] < 99) {
      if (dselect[trial] == bd) {
        bd.sel <- bd.sel + 1 / ntrial * 100
      }
      d_opt <- dselect[trial]
      if (pT.true[d_opt] > (targetT + 0.1)) {
        ov.sel <- ov.sel + 1 / ntrial * 100
      }
      if (abs(u.true[d_opt] - u.true[bd]) <= (0.05 * u.true[bd]) && d_opt <= mtd) {
        od.sel <- od.sel + 1 / ntrial * 100
      }
    }

    earlystop <- sum(dselect == 99) / ntrial * 100
    n <- rep(0, ndose)
    yE <- rep(0, ndose)
    yT <- rep(0, ndose)
    for (j in 1:ndose) {
      n[j] <- sum(as.numeric(unlist(temp$doses_given)) == j)
      yE[j] <- sum(as.numeric(unlist(temp$efficacies))
      [as.numeric(unlist(temp$doses_given)) == j])
      yT[j] <- sum(as.numeric(unlist(temp$toxicities))
      [as.numeric(unlist(temp$doses_given)) == j])
    }
    bd.pts <- bd.pts + n[bd] / ntrial / npts * 100
    dose_mask <- abs(u.true[1:mtd] - u.true[bd]) <= (0.05 * u.true[bd])
    od.pts <- od.pts + sum(n[dose_mask]) / ntrial / npts * 100
    if (n[bd] < (npts / ndose)) {
      poorall <- poorall + 1 / ntrial * 100
    }
    overdose <- overdose + sum(n[pT.true > (targetT + 0.1)]) / ntrial / npts * 100
  }

  results <- list(
    bd.sel = bd.sel, od.sel = od.sel,
    bd.pts = bd.pts, od.pts = od.pts,
    earlystop = earlystop, ntox = ntox,
    neff = neff, u.mean = u.mean,
    overdose = overdose, poorall = poorall,
    incoherent = 0, ov.sel = ov.sel
  )
  return(results)
}
