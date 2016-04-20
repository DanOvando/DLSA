#' AgeAtLength 
#' returns the age at length for a Von Bert relationship
#'
#' @param Lengths a vector of lengths
#' @param Fish a list object containing von bert parameters
#' @param Error the base standard deviation of the age at length conversion
#'
#' @return a vector of ages at length
#' @export
AgeAtLength<- function(Lengths,Fish,Error)
{
  # Error<- Fish$LengthError	
  Lengths[is.na(Lengths)]<- 0
  # Lengths<- LengthDat$Length
  AgeSD<- Error*(1+Fish$VBErrorSlope*Lengths/Fish$Linf)
  #   RawAges<- (log(1-(Lengths)/Fish$Linf)/-Fish$vbk)+Fish$t0
  RawAges<- (log(1-pmin(Lengths,Fish$Linf*.99)/Fish$Linf)/-Fish$vbk)+Fish$t0
  #   AgeWithError<- RawAges*rlnorm(length(Lengths),mean=0,sd=AgeSD)
  AgeWithError<- pmax(1,RawAges+rnorm(length(Lengths),mean=0,sd=AgeSD))
  
  return(AgeWithError)
}