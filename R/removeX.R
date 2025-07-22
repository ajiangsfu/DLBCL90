#' This is a simple function to remove "X" at the beginning of a string for a  vector of strings
#' @description  This is a simple function to remove "X" at the beginning of a string for a  vector of strings
#' @param namesIn a vector of string
#' @return a string vector
#' @author Aixiang Jiang
#' @export
removeX = function(namesIn){
  namesIn = as.character(namesIn)
  xstart = startsWith(namesIn, "X")
  namesOut = namesIn
  if(length(which(xstart)) > 0){
    namesOut[which(xstart)] = sapply(namesIn[which(xstart)], FUN = function(xx){
      substr(xx, start = 2, stop = nchar(xx))
    })
  }
  return(namesOut)
}
