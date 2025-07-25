% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lymph3Cx.R
\name{Lymph3Cx}
\alias{Lymph3Cx}
\title{A wrap-up function to get Lymph3CX calls with given model table}
\usage{
Lymph3Cx(newdat)
}
\arguments{
\item{newdat}{A new nano string data frame, samples are in columns, and genes are in rows. This is the only input parameter.
Notice that data should not be in log2 format since the log2 transformation will be done within the function}
}
\value{
A data frame with normval, PMBLval, PMBLp, PMBLcall, DLBCLval, DLBCLp, DLBCLcall, Totalcall
}
\description{
This is the function to calculate PMBL and DLBCL LPS (Linear Prediction Score) scores and empiracal Bayes' probs, 
and get DHITsig calls for a new nano data set. The required prior information is the given model table included in this package, 
which contains: 64 selected PMBL and DLBCL genes and their weights, LPS score means and sds for each of the two ends groups
and so on. The tricky is that this model table is fixed (both value and format), and the function is hard coding for the model format.
Therefore, this wrap up function loads the model table and uses it directly. 
Notice that the include.shift = T is used for the current function, if include.shift = F should be used, then should call Lymph3CX.to.DLBCL90
instead of the current wrap up function
}
\details{
Make sure the new nano data is comparable to the DLBCL90 data
}
\references{
Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A gene expression-based method
to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
A. 2003 Aug 19;100(17):9991-6.
}
\author{
Aixiang Jiang
}
\keyword{DLBCL}
\keyword{LPS}
\keyword{PMBL}
