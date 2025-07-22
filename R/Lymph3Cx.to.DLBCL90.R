#' Lymph3CX for a nano string data
#' @description This function is to get PMBL calls (DLBCL,PMBL and Unclear PMBL/GCB) and DLBCL calls (COO, cell of origin, ABC，GCB and Unclass DLBCL)
#'  for a nano string data set
#' @details Both PMBL and DLBCL classification scores are based on LPS algorithm from Wright (2003), 
#'  and their clustering probability calculation is Empiracal Bayes' based.
#' @param data a nano string data set in data frame or data matrix format, rows are for genes (symbols), and columns are for samples
#'  which should not be log transformed since log2 transformation is included in the function
#' @param model.table a known model table with fixed setting
#' @param include.shift a logic variable to indicate that if calibration should be included in the function.
#'  Its default value is TRUE, calibration is done with adding baseline shift from model table; if FALSE is chosen, however, 
#'  the calibration is reversed.
#' @return A data frame with normval, PMBLval, PMBLp, PMBLcall, DLBCLval, DLBCLp, DLBCLcall, Totalcall
#' @keywords LPS PMBL DLBCL
#' @author Stacy Hung (function), Aixiang Jiang (function description)
#' @importFrom stats dnorm
#' @references Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A gene expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.
#' @export

Lymph3CX.to.DLBCL90 <- function(data, model.table, include.shift = T)
{	#note genes in data needs to be in the same order as the example
  
  DLBCL90_CORRECTION_FACTOR = 4.4485  ### this is for PMBL
  
  # extract datasets from model table file
  paramx <- model.table[65:75,2]             # get parameters for model (last 10 lines, 2nd col of model.table.txt)
  gen.table <- model.table[1:64,]            # get model information for the genes (remainder of model.table.txt)
  
  # NB: the input data may not always have the same set of genes in the same order as the model table
  # To ensure the predictor is applied correctly:
  # 1. Remove genes from data that are not present in model table
  # 2. Remove genes from model table that are not in filtered gene table
  # 3. Keep the same order between data and model table
  
  # filter data and model table so that they have the same genes (not necessarily order)
  data <- data[row.names(data) %in% gen.table[, c("Gene.symbol")], ]
  gen.table <- gen.table[gen.table$Gene.symbol %in% row.names(data), ]
  
  # sort data and model table so they have not only same genes but also same order of genes
  data <- data[order(row.names(data)),]
  gen.table <- gen.table[order(gen.table$Gene.symbol),]
  
  dat <- logb(data,2)                        # log2 transform expression data
	
	# assign indices to model parameters
	index.norm.cut = 1
	index.DLBCL.mn = 2
	index.PMBL.mn = 3
	index.DLBCL.sd = 4
	index.PMBL.sd = 5
	index.ABC.mn = 6
	index.GCB.mn = 7
	index.ABC.sd = 8
	index.GCB.sd = 9
	index.COO.scale = 10
	index.COO.shift = 11
	
	housekeeping.genes <- unique(rbind(gen.table[grep("TRUE", gen.table$Housekeep.gen), ]))
	### all 13 house keeping genes
	
	normval <- colMeans(dat[row.names(dat) %in% housekeeping.genes$Gene.symbol, ])
	### this is used for failure test with all 13 house keeping genes, and this is mean in log2 scale
	
	PMBLval <- colSums(dat*gen.table[ ,c("pmbl.coef")])  
	     ### PMBL score calculation plus normalization with all 13 house keeping genes
	PMBLval <- PMBLval - DLBCL90_CORRECTION_FACTOR  ### calibration
	
	DLBCLval <- colSums(dat*gen.table[ ,c("DLBCL.coef")])
	     ### DLBCL score calculation plus normalization with all 5 house keeping genes (changed on Nov 17, 2020)
	DLBCLval <- DLBCLval*paramx[index.COO.scale] + paramx[index.COO.shift] ### calibration

	p1 <- dnorm(PMBLval, paramx[index.DLBCL.mn], paramx[index.DLBCL.sd])
	p2 <- dnorm(PMBLval, paramx[index.PMBL.mn], paramx[index.PMBL.sd])
	PMBLp <- p2 / (p1 + p2)
	
	PMBLcall <- 1 + (PMBLp > 0.1) + (PMBLp > 0.9)
	PMBLcall <- c("DLBCL", "Unclear", "PMBL")[PMBLcall]
	
	p1 <- dnorm(DLBCLval, paramx[index.ABC.mn], paramx[index.ABC.sd])
	p2 <- dnorm(DLBCLval, paramx[index.GCB.mn], paramx[index.GCB.sd])
	DLBCLp <- p1 / (p1 + p2)
	
	DLBCLcall <- 1 + (DLBCLp > 0.1) + (DLBCLp > 0.9)
	DLBCLcall <- c("GCB", "UNCLASS", "ABC")[DLBCLcall] ### change on 20191120 upon Barbara's request

	if(!include.shift){	
	  DLBCLval=DLBCLval-paramx[index.COO.shift]
	}	
	set=normval<paramx[index.norm.cut]
	PMBLcall[set]=DLBCLcall[set]="Low Norm"
	
	Totalcall=DLBCLcall
	set=PMBLcall=="PMBL"
	Totalcall[set]="PMBL"
	set=PMBLcall=="Unclear"
	Totalcall[set]=paste("Unclear PMBL",Totalcall[set],sep="/")
	

	data.frame(names(data),
	           normval,
	           PMBLval,
	           PMBLp,
	           PMBLcall,
	           DLBCLval,
	           DLBCLp,
	           DLBCLcall,
	           Totalcall)
	
}
