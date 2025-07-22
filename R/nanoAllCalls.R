#' A wrap-up function to get all calls: PMBL, DLBCL and DHITsig
#' @description  This is the wrap up function to get all done
#' @details This includes: nano string pre-process, PMBL + DLBCL + DHIT sig calls
#' @param nano_csv_File A nano string data file name, which can be assoicated with an absolute path or a relative path starting from your current path
#'  For example, my current path is: "/mnt/thanos_lv/ajiang/AJ2019"
#'  the file name can be: "/mnt/thanos_lv/ajiang/AJ2019/DHIT/DLBCL90 cartridges 60 to 67.csv"
#'   or: "DHIT/DLBCL90 cartridges 60 to 67.csv"
#' @param geomeanCut A QC (quality control) cutoff for house keeping genes' mean, default is 60, 
#'   which is close to 2**6 = 64. Here, 6 is the QC cutoff for PMBL/DLBCL calls that is done in log2 scale
#' @return A data frame with PMBL + DLBCL + DHIT sig calls and their related values
#' @keywords nano string 
#' @author Aixiang Jiang
#' @export
nanoAllCalls = function(nano_csv_File, geomeanCut = 60){
   dat = nanoProcess(nano_csv_File)
   dlbcl = Lymph3Cx(dat$COOdat)
   dhits = DHITsig(dat$DHITdat)
   ### add QC step:  QC step (already in place for COO/PMBCL where the call is FAIL or LOW NORM instead of a definitive call for DHITsig 
   ### if the geomean is less than 64 (same threshold as for the rest of the algorithm)
   tmp = which(dat$house_geomean < geomeanCut)
   dhits[tmp, "DHITsig_class"] = "Low Norm"
   
   allres = cbind(dlbcl, dhits[rownames(dlbcl),])
   
   ### make changes on 20190325
   ### change format first
   tmp = grep("call", colnames(allres))
   ### add DHITsig_class in with tmp, for each of these columns, do the following
   
   allres[,c(tmp,11)] = apply(allres[,c(tmp,11)], 2, function(xx){
     xx=as.character(xx)
     xx = gsub("Low Norm", "Fail",xx)
   })
   
   ### for Totalcall, do the following
   #   if PMBLcall == PMBL:
   #   totalcall = PMBLcall
   # else:
   #   totalcall = DLBCLcall (edited) 
   
   allres$Totalcall = ifelse(allres$PMBLcall == "PMBL", "PMBL",allres$DLBCLcall)
   
   return(allres)
}
