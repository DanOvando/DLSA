#' Check Sample Size
#'
#' Checks the sample size of the length data to make sure you've got enough
#'
#' @param Data
#'
#' @return
#' @export
CheckLengthSampleSize<- function(Data)
{
  #   Data<- LengthData
  SampleSize<- Data %>%
    group_by(Year) %>%
      summarize(SampleSize=length(Length))

  SufficientSamples<- SampleSize$SampleSize>MinSampleSize

  if (sum(SufficientSamples)>0)
  {
    AvailableYears<- SampleSize$Year[SufficientSamples]

    Data<- Data[Data$Year %in% AvailableYears,]
  }

  return(list(YearsWithEnoughData=sum(SufficientSamples),ParedData=Data))
}
