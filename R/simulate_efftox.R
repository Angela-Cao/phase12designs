#' Simulate operating characteristics using EffTox
#'
#' This function runs simulations of the EffTox design by
#' evaluating operating characteristics over a range of cohort sizes. For each
#' dose level within the user-specified range, it performs multiple trials and saves the results to a corresponding file.
#'
#' @param ndose Integer. Number of dose levels. (**Required**)
#' @param ssizerange Integer vector. Range of number of cohorts to simulate. (**Required**)
#' @param target_t Numeric. Target toxicity probability. (**Required**)
#' @param lower_e Numeric. Minimum acceptable efficacy probability. (**Required**)
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
#' simulate_efftox(
#'   ndose = 5,
#'   ssizerange = 1:2,
#'   target_t = 0.3,
#'   lower_e = 0.4,
#'   ntrial = 10,
#'   prob = prob,
#' )

#' @export
simulate_efftox <- function(ndose, ssizerange, target_t, lower_e,
                            startdose = 1, ntrial = 10000,
                            utilitytype = 1, prob = NULL, save_dir = tempdir(),
                            save_folder = "efftox_simulations",
                            save_file = "efftox_simulation.csv") {
  full_save_root <- file.path(save_dir, save_folder)
  dir.create(full_save_root, recursive = TRUE, showWarnings = FALSE)

  for (iii in 1:ndose) {
    OBD <- iii

    xmat <- NULL
    for (utype in c(1, 2)) {
      for (rtype in c(1)) {
        for (i in ssizerange) {
          xmat <- rbind(xmat, c(i, utype, rtype))
        }
      }
    }

    outputmat <- NULL
    ooo <- NULL

    for (kk in 1:nrow(xmat)) {
      i <- xmat[kk, 1]
      utype <- xmat[kk, 2]
      rtype <- xmat[kk, 3]

      oc <- oc_efftox(
        ndose = ndose, target_t = target_t, lower_e = lower_e,
        ncohort = i, startdose = startdose, ntrial = ntrial,
        utilitytype = utype, prob = prob,OBD = OBD
      )
      # print(i)
      result_row <- c(
        i, utype, rtype,
        oc$bd.sel, oc$od.sel, oc$bd.pts, oc$od.pts,
        oc$earlystop, oc$overdose, oc$poorall, oc$ov.sel
      )

      outputmat <- rbind(outputmat, result_row)
    }

    outputmat <- as.data.frame(outputmat)
    colnames(outputmat) <- c(
      "ncohort", "utype", "rtype",
      "bd.sel", "od.sel", "bd.pts", "od.pts",
      "earlystop", "overdose", "poorall", "ov.sel"
    )

    outputmat$design <- "efftox"

    subfolder_path <- file.path(full_save_root, as.character(iii))
    dir.create(subfolder_path, showWarnings = FALSE)
    file_path <- file.path(subfolder_path, save_file)
    write.csv(outputmat, file_path, row.names = FALSE)
  }
}
