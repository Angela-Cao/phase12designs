#' Simulate operating characteristics using PRINTE
#'
#' This function runs simulations of the PRINTE design by
#' evaluating operating characteristics over a range of cohort sizes. For each
#' dose level within the user-specified range, it performs multiple trials and saves the results to a corresponding file.
#'
#' @param ndose Integer. Number of dose levels. (**Required**)
#' @param ssizerange Integer vector. Range of number of cohorts to simulate. (**Required**)
#' @param target_t Numeric. Target toxicity probability. (**Required**)
#' @param target_e Numeric. Target efficacy probability. (**Required**)
#' @param lower_e Numeric. Minimum acceptable efficacy probability. (**Required**)
#' @param cohortsize Integer. Size of a cohort. (Default is `3`)
#' @param startdose Integer. Starting dose level. (Default is `1`)
#' @param eps1 Numerical. Width of the subrectangle.
#' @param eps2 Numerical. Width of the subreactangle.
#' @param psafe Numeric. Early stopping cutoff for toxicity. (Default is `0.95`)
#' @param pfutility Numeric. Early stopping cutoff for efficacy. (Default is `0.90`)
#' @param ntrial  Integer. Number of random trial replications. (Default is `10000`)
#' @param utilitytype Integer. Type of utility structure. (Default is `1`)
#'   - If set to `1`: Use preset weights (w11 = 0.6, w00 = 0.4)
#'   - If set to `2`: Use (w11 = 1, w00 = 0)
#'   - Other: Use user-specified values from `u1` and `u2`.
#' @param u1 Numeric. Utility parameter w_11. (0-100)
#' @param u2 Numeric. Utility parameter w_00. (0-100)
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
#' @param save_dir Directory to save output folders. Default is `tempdir()`.
#' @param save_folder Folder name. (Default is "boin12_simulations")
#' @param save_file File name. (Default is "boin12_simulation.csv")
#' @return No return value, called for side effects
#' @examples
#' prob <- list(
#'   pE = c(0.4, 0.5, 0.6, 0.6, 0.6),
#'   pT = c(0.1, 0.2, 0.3, 0.4, 0.4),
#'   obd = 3,
#'   mtd = 2
#' )
#' simulate_pite(
#'   ndose = 5,
#'   ssizerange = 1:2,
#'   target_t = 0.3,
#'   target_e = 0.5,
#'   lower_e = 0.4,
#'   ntrial = 10,
#'   prob = prob,
#' )

#' @export
simulate_pite <- function(ndose, ssizerange,
                          target_t, target_e, lower_e,
                          cohortsize = 3, startdose = 1,
                          eps1 = 0.05, eps2 = 0.05,
                          psafe = 0.95, pfutility = 0.9,
                          ntrial = 10000, utilitytype = 1,
                          u1, u2, prob = NULL,
                          save_dir = tempdir(), save_folder = "pite_simulations",
                          save_file = "pite_simulation.csv") {
  full_save_root <- file.path(save_dir, save_folder)
  dir.create(full_save_root, recursive = TRUE, showWarnings = FALSE)

  for (iii in 1:ndose) {
    OBD <- iii
    outputmat <- NULL
    for (utype in c(1, 2)) {
      for (rtype in c(1)) {
        for (i in ssizerange) {
          oc <- oc_pite(
            ndose = ndose, target_t = target_t, target_e = target_e,
            lower_e = lower_e, ncohort = i, cohortsize = cohortsize,
            startdose = startdose, eps1 = eps1, eps2 = eps2,
            psafe = psafe, pfutility = pfutility,
            ntrial = ntrial,OBD = OBD,
            utilitytype = utype, u1 = u1, u2 = u2, prob = prob
          )
          outputmat <- rbind(outputmat, c(i, utype, rtype, c(oc$bd.sel, oc$od.sel, oc$bd.pts, oc$od.pts, oc$earlystop, oc$overdose, oc$poorall, oc$ov.sel)))
        }
      }
    }
    cname <- c("ncohort", "utype", "rtype", "bd.sel", "od.sel", "bd.pts", "od.pts", "earlystop", "overdose", "poorall", "ov.sel")
    colnames(outputmat) <- cname
    cname <- c(cname, "design")
    outputmat <- cbind(outputmat, rep("printe", nrow(outputmat)))
    colnames(outputmat) <- cname

    subfolder_path <- file.path(full_save_root, as.character(iii))
    dir.create(subfolder_path, showWarnings = FALSE)
    file_path <- file.path(subfolder_path, save_file)
    write.csv(outputmat, file_path, row.names = FALSE)
  }
}
