% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/source_GESTIA.R
\name{perturb_pwPos}
\alias{perturb_pwPos}
\title{Permutation of the inter-pathway interactions}
\usage{
perturb_pwPos(g.ab, gss1, gss2)
}
\arguments{
\item{g.ab}{The network of the two pathways}

\item{gss1}{The gene symbols of pathway A's genes}

\item{gss2}{The gene symbols of pathway B's genes}
}
\value{
An igraph object of the permutated network of A and B
}
\description{
Randomly assign the edges between two pathways, maintains the number of edges in/out the pathways, remains the structure of each pathway
}
