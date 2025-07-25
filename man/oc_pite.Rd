% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oc_pite.R
\name{oc_pite}
\alias{oc_pite}
\title{Compute Operating Characteristics using PRINTE}
\usage{
oc_pite(
  ndose,
  target_t,
  target_e,
  lower_e,
  ncohort = 10,
  cohortsize = 3,
  startdose = 1,
  eps1 = 0.05,
  eps2 = 0.05,
  psafe = 0.95,
  pfutility = 0.95,
  ntrial = 10000,
  utilitytype = 1,
  u1,
  u2,
  prob = NULL
)
}
\arguments{
\item{ndose}{Integer. Number of dose levels. (\strong{Required})}

\item{target_t}{Numeric. Target toxicity probability. (\strong{Required})}

\item{target_e}{Numeric. Target efficacy probability. (\strong{Required})}

\item{lower_e}{Numeric. Minimum acceptable efficacy probability. (\strong{Required})}

\item{ncohort}{Integer. Number of cohorts. (Default is \code{10})}

\item{cohortsize}{Integer. Size of a cohort. (Default is \code{3})}

\item{startdose}{Integer. Starting dose level. (Default is \code{1})}

\item{eps1}{Numerical. Width of the subrectangle. (Default is '0.05')}

\item{eps2}{Numerical. Width of the subreactangle. (Default is '0.05')}

\item{psafe}{Numeric. Early stopping cutoff for toxicity. (Default is \code{0.95})}

\item{pfutility}{Numeric. Early stopping cutoff for efficacy. (Default is \code{0.95})}

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
}
\value{
A list containing operating characteristics such as:
\describe{
\item{bd.sel}{OBD selection percentage}
\item{od.sel}{Favorable dose selection percentage}
\item{bd.pts}{Average percentage of patients at the OBD }
\item{od.pts}{Average percentage of patients at the favorable doses}
\item{earlystop}{Percentage of early stopped trials}
\item{overdose}{Overdose patients percentage }
\item{poorall}{Poor allocation percentage}
\item{ov.sel}{Overdose selection percentage}
}
}
\description{
\code{oc_pite()} uses the PRINTE design to compute operating charateristics of a user-specificed trial scenario.
This design maps toxicity and efficacy intervals onto a decision table, forming 16 equal-area regions.
}
