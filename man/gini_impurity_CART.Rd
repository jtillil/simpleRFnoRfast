% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{gini_impurity_CART}
\alias{gini_impurity_CART}
\title{Gini impurity of candidate linear splits in CART algorithm.}
\usage{
gini_impurity_CART(dat, label, pure_coefs, value, candidates, varID, gamma)
}
\arguments{
\item{dat}{Data frame of predictors.}

\item{label}{Vector of labels.}

\item{par}{Vector of coefficients of split.}
}
\value{
Gini impurities.
}
\description{
Compute Gini impurity of candidate linear splits in CART algorithm.
}
\author{
Johannes Tillil
}
