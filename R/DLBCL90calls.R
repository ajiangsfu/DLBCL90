#' A wrap-up function to get all calls: PMBL, DLBCL and DHITsig with two different data formats
#' @description  This is the wrap up function to get all done and for input data with either 1 and more heading rows.
#' @details This includes: nano string data re-formatring when necessary, nano string pre-process, PMBL + DLBCL + DHIT sig calls, and combination
#' @param DLBCL90_File A nano string data file name, which can be assoicated with an absolute path or a relative path starting from your current path,
#'   it can be either raw data as a zip file name, or a combined data file in csv format.
#' @param nHeading A integer of number of heading rows, default is 2 since most recently files are with two rows of heading info, however, 
#'   I only check if this is 1 or not
#' @param outfileName An output file name, which can be assoicated with an absolute path or a relative path starting from your current path
#' @param geomeanCut A QC (quality control) cutoff for house keeping genes' mean, default is 60, 
#'   which is close to 2**6 = 64. Here, 6 is the QC cutoff for PMBL/DLBCL calls that is done in log2 scale
#' @return A data frame with PMBL + DLBCL + DHIT sig calls and their related values. At the same time, this function writes result table 
#'   into the same path of input data, also writes one row heading data table as well heading info table if input file contains two heading rows.
#' @keywords DLBCL90 nano-string 
#' @importFrom utils read.csv write.csv
#' @author Aixiang Jiang
#' @export

DLBCL90calls = function(DLBCL90_File, outfileName, nHeading = 2, geomeanCut = 60){
  
  #### first of all, check if the input file is a zip file name or a csv file name
  #### if it is a zip file name, use nanostringr::read_rcc function to deal with zip file
  #### if the giving file containing any path info, I should extract it from the zip file name, and put it inot read_rcc function
  #### 
  
  # csvtype = ifelse(endsWith(DLBCL90_File, "CSV") | endsWith(DLBCL90_File, "csv"), 1, 0)
  
  ziptype = ifelse(endsWith(DLBCL90_File, "ZIP") | endsWith(DLBCL90_File, "zip"), 1, 0)
  
  ### change on 20200806
  fullfiles = NA
  
  if(ziptype == 1){
    ### get the path, break from the last "/"
    tmp = strsplit(DLBCL90_File, split = "/")[[1]]
    tn = length(tmp)
    pathName = "."
    if(tn > 1){
      tmp = paste("/",tmp[tn], sep="")
      pathName = gsub(tmp, "", DLBCL90_File)
    }
    
    ### unzip before doing anything else in the same folder
    utils::unzip(zipfile = DLBCL90_File, exdir = pathName)
    
    ### change on 20200806
    fullfiles = list.files(path = pathName, pattern = "\\.RCC$", full.names = TRUE, 
                           ignore.case = TRUE)
    
    filesIn = nanostringr::read_rcc(path = pathName)
    csvfile = data.frame(filesIn$raw)
    rownames(csvfile) = csvfile$Name
    csvfile = csvfile[,-c(1:3)]
    
    ### to avoid log(0) is inf, for all 0, change to 1
    tmp = which(csvfile < 1, arr.ind = TRUE)
    if(dim(tmp)[1] > 0) {csvfile[tmp] = 1}
    
    nn = nchar(DLBCL90_File)
    fileout = substr(DLBCL90_File, start = 1, stop = nn-4)
    fileout = paste(fileout, ".csv", sep="")
    colnames(csvfile) = removeX(colnames(csvfile))
    write.csv(csvfile, fileout)
    DLBCL90_File = fileout
    nHeading = 3
  }
  
  allout = NA
  
  if (nHeading == 1) {
    allout = nanoAllCalls(nano_csv_File = DLBCL90_File, geomeanCut = geomeanCut)
    colnames(allout)[1] = "sampleName"
    allout$sampleName = removeX(allout$sampleName)
    
  } else if (nHeading == 2){
    workdat = conv2NamesFile(fileIn = DLBCL90_File)
    calls = nanoAllCalls(workdat[1])
    info = read.csv(
      workdat[2],
      header = T,
      row.names = 1,
      stringsAsFactors = F
    )
    allout = cbind(info, calls)
    allout = allout[, -2]
    allout$longName = rownames(allout)
    allout = allout[, c(14, 1:13)]
    allout$longName = removeX(allout$longName)
  } else {  ### add this choice on 20200806
    allout = nanoAllCalls(nano_csv_File = DLBCL90_File, geomeanCut = geomeanCut)
    colnames(allout)[1] = "sampleName"
    allout$sampleName = removeX(allout$sampleName)
    
    ### modify fullfiles
    fullfiles = sapply(fullfiles, function(xx){
      tmp = strsplit(xx, split = "/")[[1]]
      nn = length(tmp)
      return(tmp[nn])
    })
    
    allout$longName = fullfiles
    allout = allout[,c(14,1:13)]
    colnames(allout)[2] = "shortName"
  }
  
  
  rownames(allout) = removeX(rownames(allout))
  
  ### make sure the outfileName is ended with ".csv"
  nn = nchar(outfileName)
  ends = substr(outfileName, start = nn-3, stop = nn)
  if(ends != ".csv"){
    outfileName = paste(outfileName, ".csv", sep="")
  }
  
  write.csv(allout, outfileName, row.names = F)
  
  return(allout)
}

