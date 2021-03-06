% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppR6.R, R/es20.R
\name{ES20_Environment}
\alias{ES20_Environment}
\alias{ES20_make_environment}
\alias{ES20_fixed_environment}
\title{ES20_Environment object}
\usage{
ES20_Environment(disturbance_mean_interval, seed_rain, k_I, control)

ES20_make_environment(p)

ES20_fixed_environment(e = 1, p = ES20_Parameters(), height_max = 150)
}
\arguments{
\item{disturbance_mean_interval}{mean disturbance interval in years}

\item{seed_rain}{vector of seed rains}

\item{k_I}{light extinction coefficient}

\item{control}{A Control object}

\item{p}{A Parameters object}

\item{height_max}{= 150.0 maximum possible height in environment}

\item{e=1.0}{Value of environment}
}
\description{
ES20_Environment object

Construct a fixed environment for ES20 strategy
}
