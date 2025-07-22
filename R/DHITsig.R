#' A wrap-up function to get DHITsig calls 
#' @description This is the function to calculate DHITsig LPS (Linear Predictor Score) scores and empiracal Bayes' probs, 
#' and get DHITsig calls for a new nano data set. The required prior information is: 
#' 30 selected DHITsig genes and their weights, LPS score means and sds for DHITsig POS and NEG groups 
#' @details Make sure the new nano data is comparable to the DLBCL90 data, if not, calibration is required before calling this function
#' @param newdat A new nano string data frame, samples are in columns, and genes are in rows. 
#'  Notice that data are already pre-processed, normalized with house keeping genes, and in log2 format
#' @param classProbCut A number for probability cutoff, which should be within (0.5, 1), the default value is 0.8
#' @return A data frame with LPS score, Empirical Bayesian probabilites for two groups and classification
#' @keywords DHITsig, LPS, nano string
#' @author Aixiang Jiang
#' @references 
#' Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A gene expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.
#' 
#' Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S, Rushton C, Tang J, Thomas N, Slack GW, Farinha P, 
#' Takata K, Miyata-Takata T, Craig J, Mottok A, Meissner B, Saberi S, Bashashati A, Villa D, Savage KJ, Sehn LH, Kridel R,
#' Mungall AJ, Marra MA, Shah SP, Steidl C, Connors JM, Gascoyne RD, Morin RD, Scott DW. 
#' Double-Hit Trait Expression Signature Defines a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse Large B-Cell
#' Lymphoma. J Clin Oncol. 2018 Dec 3:JCO1801583. doi: 10.1200/JCO.18.01583.
#' @export

DHITsig = function(newdat,  classProbCut = 0.8){  ## notice that the cutoff is 0.8 for DLC90 DHITsig calls
  ## assuming that the package is loaded, and we can get the following two data directly
  load(system.file("extdata", "DHITsig30nano.rda", package = "DLBCL90"))
  load(system.file("extdata", "DHITsigLPSpars.rda", package = "DLBCL90"))
  outs = LPSClassWithPriors(newdat = newdat, weights = DHITsig30nano,testGroup = "POS", refGroup = "NEG", 
                                LPSMeanSds = c(DHITsigLPSpars[,1], DHITsigLPSpars[,2]), classProbCut = classProbCut)
  colnames(outs) = gsub("LPS","DHITsig", colnames(outs))
  colnames(outs) = gsub("test", "POS", colnames(outs))
  colnames(outs) = gsub("ref", "NEG", colnames(outs))
  return(outs)
}














