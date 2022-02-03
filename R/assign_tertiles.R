#' assign_tertiles
#' Divides continuous gene expression levels into 3 tertiles and label them "_high", "_medium" and "_low" concatenating with the input id text.
#' @param numeric_vector1 any continuous numeric vector
#' @param id1 A character of length 1. Used to describe the tertiles
#' @return character vector of the same length with input numeric vector. The values are replaced by 3 levels of tertiles.
#'
#' @author Tolga Turan, \email{tolga.turan@abbvie.com}
#' @references \url{http://github.com/tolgaturan-github/IBRIDGE}
#' @examples
#' tertiles1<-assign_tertiles(numeric_vector1, "inflamed")




assign_tertiles<-function(numeric_vector1, id1){
	n<-cut(numeric_vector1,quantile(numeric_vector1, 0:3/3,na.rm=TRUE), include.lowest=TRUE)
	levels(n)<-c(paste0(id1,"_low"), paste0(id1,"_medium"), paste0(id1,"_high"))
	n}
