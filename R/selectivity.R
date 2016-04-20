#' Title
#'
#' @param length
#' @param shape
#' @param form
#'
#' @return return selectivity at length ogive
#' @export

selectivity <- function(length,shape,form = 'logistic'){

  if (form == 'logistic'){
  selectivity_ogive <- 1.0/(1+ exp(-log(19)*(length - shape$SL50)/(shape$SL95 - shape$SL50)))
  }

  return(selectivity_ogive)
}
