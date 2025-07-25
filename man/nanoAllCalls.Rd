% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nanoAllCalls.R
\name{nanoAllCalls}
\alias{nanoAllCalls}
\title{A wrap-up function to get all calls: PMBL, DLBCL and DHITsig}
\usage{
nanoAllCalls(nano_csv_File, geomeanCut = 60)
}
\arguments{
\item{nano_csv_File}{A nano string data file name, which can be assoicated with an absolute path or a relative path starting from your current path
For example, my current path is: "/mnt/thanos_lv/ajiang/AJ2019"
the file name can be: "/mnt/thanos_lv/ajiang/AJ2019/DHIT/DLBCL90 cartridges 60 to 67.csv"
 or: "DHIT/DLBCL90 cartridges 60 to 67.csv"}

\item{geomeanCut}{A QC (quality control) cutoff for house keeping genes' mean, default is 60, 
which is close to 2**6 = 64. Here, 6 is the QC cutoff for PMBL/DLBCL calls that is done in log2 scale}
}
\value{
A data frame with PMBL + DLBCL + DHIT sig calls and their related values
}
\description{
This is the wrap up function to get all done
}
\details{
This includes: nano string pre-process, PMBL + DLBCL + DHIT sig calls
}
\author{
Aixiang Jiang
}
\keyword{nano}
\keyword{string}
