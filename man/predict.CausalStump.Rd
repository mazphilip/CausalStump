% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.CausalStump.R
\name{predict.CausalStump}
\alias{predict.CausalStump}
\title{Predict the response surfaces using a Gaussian process or Student-t process fit.}
\usage{
\method{predict}{CausalStump}(cs_object, X, z, pscore, nsampling)
}
\arguments{
\item{cs_object}{A object from the CausalStump function}

\item{X}{A data.frame with new data. If not presented, using the training data.}

\item{z}{A vector or scalar with new treatment data. If it is a scalar, it predicts using the same value for all observations. If missing, it uses the training data.}

\item{pscore}{A vector with the propensity score. Throws an error if use is different to the CausalStump fit.}

\item{nsampling}{(optional) A number that overwrites the number of samples used for the TP-kernel sampling.}
}
\value{
Returns the MAP and the 95 percent credible interval of the fitted process as a list.
}
\description{
Predict the response surfaces using a Gaussian process or Student-t process fit.
}
