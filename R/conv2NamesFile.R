#' Processing data with two or more heading rows
#' @description  This is to separate input data containing two or more heading rows into regular data format table and one info table
#' @details This function reads in full path csv file with two or more rows of heading info, then write out a usual one row heading data table and an 
#'  info table contain all heading info in the same path of the input file.
#' @param fileIn An input nano string data file name, which can be assoicated with an absolute path or a relative path starting from your current path.
#'   The file itself contains two or more rows of heading, which is not usually data format for R, but we need it for nano string DLBCL90 to make
#'   sure we have enough information
#' @return A string vector containing two items, the 1st one is the regular data format name, and the 2nd one is the file name with extra sample info
#' @keywords nano string 
#' @author Aixiang Jiang
#' @export


conv2NamesFile = function(fileIn){
  
  rawdat = read.csv(fileIn, header=T, row.names = 1, stringsAsFactors = F)

  datfile = rawdat
  namefile = data.frame(t(rawdat), stringsAsFactors = F)
  tmp = namefile
  
  t1 = which(rownames(datfile) == "")  ### this can deal with multiple extra headings, whether it is one or more extra rows
  if(length(t1) > 0) {
    datfile = datfile[-t1,]
    namefile = namefile[,t1]
  }

  namefile = data.frame(namefile,stringsAsFactors = F)
  
  ### remove NA
  namefile = na.omit(namefile)
  tmp = na.omit(tmp)
  
  rownames(namefile) = rownames(tmp)
  colnames(namefile)[1] = "shortName"  ### only re-define the 1st extra description, although in theory, multiple extra descriptions should work as well
  
  
  ### now, need to remove Annotation and any other empty columns
  datfile = datfile[,rownames(namefile)]
  
  datfile = apply(datfile, c(1,2), as.numeric)
  datfile = data.frame(datfile, stringsAsFactors = F)
  
  fileIn = gsub(".csv", "", fileIn)
  datOut = paste(fileIn, "_datPart.csv", sep="")
  nameOut = paste(fileIn, "_NamePart.csv", sep="")
  
  write.csv(datfile, datOut)
  write.csv(namefile, nameOut)
  
  return(c(datOut, nameOut))
  
}

