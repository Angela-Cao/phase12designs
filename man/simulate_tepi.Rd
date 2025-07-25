% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_tepi.R
\name{simulate_tepi}
\alias{simulate_tepi}
\title{Simulate Operating Characteristics using TEPI}
\usage{
simulate_tepi(
  ndose,
  ssizerange,
  target_t,
  lower_e,
  cohortsize = 3,
  startdose = 1,
  effint_l = c(0, lower_e, lower_e + 0.2, lower_e + 0.4),
  effint_u = c(lower_e, lower_e + 0.2, lower_e + 0.4, 1),
  toxint_l = c(0, 0.15, target_t, target_t + 0.05),
  toxint_u = c(0.15, target_t, target_t + 0.05, 1),
  psafe = 0.95,
  pfutility = 0.9,
  ntrial = 10000,
  utilitytype = 1,
  u1,
  u2,
  prob = NULL,
  save_dir = ".",
  save_folder = "tepi_simulations",
  save_file = "tepi2_simulation.csv"
)
}
\arguments{
\item{ndose}{Integer. Number of dose levels. (\strong{Required})}

\item{ssizerange}{Integer vector. Range of number of cohorts to simulate. (\strong{Required})}

\item{target_t}{Numeric. Target toxicity probability. (\strong{Required})}

\item{lower_e}{Numeric. Minimum acceptable efficacy probability. (\strong{Required})}

\item{cohortsize}{Integer. Size of a cohort. (Default is \code{3})}

\item{startdose}{Integer. Starting dose level. (Default is \code{1})}

\item{effint_l}{Lower efficacy bounds for dose assignment decision table. (Default is \code{c(0,lower_e,lower_e+0.2,lower_e+0.4)})}

\item{effint_u}{Lower efficacy bounds for dose assignment decision table. (Default is \code{c(lower_e,lower_e+0.2,lower_e+0.4,1)})}

\item{toxint_l}{Lower toxicity bounds for dose assignment decision table. (Default is \code{c(0,0.15,target_t,target_t+0.05)})}

\item{toxint_u}{Lower toxicity bounds for dose assignment decision table. (Default is \code{c(0.15,target_t,target_t+0.05,1)})}

\item{psafe}{Numeric. Early stopping cutoff for toxicity. (Default is \code{0.95})}

\item{pfutility}{Numeric. Early stopping cutoff for efficacy. (Default is \code{0.90})}

\item{ntrial}{Integer. Number of random trial replications. (Default is \code{10000})}

\item{utilitytype}{Integer. Type of utility structure. (Default is \code{1})
\itemize{
\item If set to \code{1}: Use preset weights (w11 = 0.6, w00 = 0.4)
\item If set to \code{2}: Use (w11 = 1, w00 = 0)
\item Other: Use user-specified values from \code{u1} and \code{u2}.
}}

\item{u1}{Numeric. Utility parameter w_11. (0-100)}

\item{u2}{Numeric. Utility parameter w_00. (0-100)}

\item{prob}{Fixed probability vectors. If not specified, a random scenario is used by default.
Use this parameter to provide fixed probability vectors as a list with the following named elements:
\itemize{
\item \code{pE}: Numeric vector of efficacy probabilities for each dose level.
\item \code{pT}: Numeric vector of toxicity probabilities for each dose level.
\item \code{obd}: Integer indicating the index of the true Optimal Biological Dose (OBD).
\item \code{mtd}: Integer indicating the index of the true Maximum Tolerated Dose (MTD).
}

For example:

\if{html}{\out{<div class="sourceCode r">}}\preformatted{prob <- list(
  pE = c(0.4, 0.5, 0.6, 0.6, 0.6),
  pT = c(0.1, 0.2, 0.3, 0.4, 0.4),
  obd = 3,
  mtd = 2
)
}\if{html}{\out{</div>}}}

\item{save_dir}{Directory to save output folders. Default is (\code{"."}).}

\item{save_folder}{Folder name. (Default is "boin12_simulations")}

\item{save_file}{File name. (Default is "boin12_simulation.csv")}
}
\value{
Results are saved as CSV files organized by OBD within folders.
}
\description{
This function runs simulations of the TEPI design by
evaluating operating characteristics over a range of cohort sizes. For each
dose level within the user-specified range, it performs multiple trials and saves the results to a corresponding file.
}
