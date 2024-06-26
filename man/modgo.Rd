% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modgo.R
\name{modgo}
\alias{modgo}
\title{MOck Data GeneratiOn}
\usage{
modgo(
  data,
  ties_method = "max",
  variables = colnames(data),
  bin_variables = NULL,
  categ_variables = NULL,
  count_variables = NULL,
  n_samples = nrow(data),
  sigma = NULL,
  nrep = 100,
  noise_mu = FALSE,
  pertr_vec = NULL,
  change_cov = NULL,
  change_amount = 0,
  seed = 1,
  thresh_var = NULL,
  thresh_force = FALSE,
  var_prop = NULL,
  var_infl = NULL,
  infl_cov_stable = FALSE,
  tol = 1e-06,
  stop_sim = FALSE,
  new_mean_sd = NULL,
  multi_sugg_prop = NULL,
  generalized_mode = FALSE,
  generalized_mode_model = NULL,
  generalized_mode_lmbds = NULL
)
}
\arguments{
\item{data}{A data frame containing the data whose characteristics are to be
mimicked during the data simulation.}

\item{ties_method}{Method used to deal with ties
during rank transformation. Allowed input: "max","average" or "min". This
parameter is passed by \code{\link[modgo]{rbi_normal_transform}} to the
parameter \code{ties.method} of \code{\link[base]{rank}}.}

\item{variables}{A character vector indicating the columns in \code{data} 
to be used. Default: \code{colnames(data)}.}

\item{bin_variables}{A character vector listing those entries in 
\code{variables} to be treated as binary variables.}

\item{categ_variables}{A character vector listing those entries in 
\code{variables} to be treated as ordinal categorical variables, with 
more than two categories. See Details.}

\item{count_variables}{A character vector listing those entries 
\code{categ_variables}  to be treated as count variables. Relevant only when
\code{generalized_mode = TRUE}.}

\item{n_samples}{Number of rows of each simulated dataset. Default is
the number of rows of \code{data}.}

\item{sigma}{A covariance matrix of NxN (N= number of variables)
provided by the user to bypass the covariance matrix calculations}

\item{nrep}{Number of simulated datasets to be generated.}

\item{noise_mu}{Logical. Should noise be added to the
mean vector of the multivariate normal distribution used to draw the 
simulated values? Default: FALSE.}

\item{pertr_vec}{A named vector. Vector's names are the continuous variables
that the user want to perturb. Variance of simulated dataset mimic original
data's variance.}

\item{change_cov}{Change the covariance of a specific pair of variables.}

\item{change_amount}{the amount of change in  the covariance
of a specific pair of variables.}

\item{seed}{A numeric value specifying the random seed. If \code{seed = NA},
no random seed is set.}

\item{thresh_var}{A data frame that contains the thresholds(left and right)
of specified variables
(1st column: variable names, 2nd column: Left thresholds,
3rd column: Right thresholds)}

\item{thresh_force}{A logical value indicating if you want to force threshold
in case the proportion of samples that can surpass the threshold are less
than 10\%}

\item{var_prop}{A named vector that provides a  proportion of
value=1 for a specific binary variable (=name of the vector) that will be
the proportion of this value in the simulated datasets.[this may increase
execution time drastically]}

\item{var_infl}{A named vector. Vector's names are the continuous variables
that the user want to perturb and increase their variance}

\item{infl_cov_stable}{Logical value. If TRUE,perturbation is applied to
original dataset and simulations values mimic the perturbed original 
dataset. Covariance matrix used for simulation = original data's correlations.
If FALSE, perturbation is applied to the simulated datasets.}

\item{tol}{A numeric value that set up
tolerance(relative to largest variance) for numerical lack of
positive-definiteness in Sigma}

\item{stop_sim}{A logical value indicating if the analysis should
stop before simulation and produce only the correlation matrix}

\item{new_mean_sd}{A matrix that contains two columns named
"Mean" and "SD" that the user specifies desired Means and Standard Deviations
in the simulated datasets for specific continues variables. The variables
must be declared as ROWNAMES in the matrix.}

\item{multi_sugg_prop}{A named vector that provides a  proportion of
value=1 for specific binary variables (=name of the vector) that will be
the close to the proportion of this value in the simulated datasets.}

\item{generalized_mode}{A logical value indicating if generalized lambda/Poisson
distributions or set up thresholds will be used to generate the simulated values}

\item{generalized_mode_model}{A matrix that contains two columns named 
"Variable" and "Model". This matrix can be used only if a 
generalized_mode_model argument is provided. It specifies what model should 
be used for each Variable. Model values should be "rmfmkl", "rprs", "star" 
or a combination of them, e.g. "rmfmkl-rprs" or "star-star", in case the 
user wants a bimodal simulation.
The user can select Generalised Poisson model for Poisson variables,
but this model cannot be included in bimodal simulation}

\item{generalized_mode_lmbds}{A matrix that contains lambdas values for each 
of the variables of the dataset to be used for either Generalized Lambda 
Distribution Generalized Poisson Distribution or setting up thresholds}
}
\value{
A list with the following components:
\item{simulated_data}{A list of data frames containing the simulated data.}
\item{original_data}{A data frame with the input data.}
\item{correlations}{A list of correlation matrices. The ith element is the
correlation matrix for the ith simulated dataset. The \code{(repn + 1)}the
(last) element of the list is the average of the correlation matrices.}
\item{bin_variables}{A character vector listing the binary variables}
\item{categ_variables }{A character vector listing the ordinal
categorical variables}
\item{covariance_matrix}{Covariance matrix used when generating observations
from a multivariate normal distribution.}
\item{seed}{Random seed used.}
\item{samples_produced}{Number of rows of each simulated dataset.}
\item{sim_dataset_number}{Number of simulated datasets produced.}
}
\description{
Creates synthetic dataset based on real data by means of the
rank-based inverse normal transformation. Data with perturbed
characteristics can be generated.
}
\details{
Simulated data is generated based on available data. The algorithm 
used is based on the ranked-based inverse normal transformation 
(Koliopanos et al. (2023)) and attempts to mimic the characteristics of the 
original data. 

All variables passed to \code{modgo} should be of class 
double or integer. This includes the variables passed to the parameter
\code{categ_variables}. The character vector \code{variables}, indicating
the variables in \code{data} to be used in the simulation, should 
contain at least two variables. The variables in \code{variables} not present
in \code{bin_variables} nor \code{categ_variables} will be treated as 
continuous variables.
}
\examples{
data("Cleveland",package="modgo")
test_modgo <- modgo(data = Cleveland,
     bin_variables = c("CAD","HighFastBloodSugar","Sex","ExInducedAngina"),
     categ_variables =c("Chestpaintype"))
}
\references{
Koliopanos, G., Ojeda, F. and Ziegler A. (2023).
A simple-to-use R package for mimicking study data by simulations.
\emph{Methods Inf Med}, 62(03/04), 119-129.
}
\author{
Francisco M. Ojeda, George Koliopanos
}
\keyword{data}
\keyword{generation}
\keyword{mock}
