#' AddMissingFish adds in zero observations for missing fish
#'
#' @param Data
#'
#' @return data with missing observations added in
#' @export
AddMissingFish<- function(Data) {

  #   Data<- GFD
  LifeStart<- which(colnames(Data)=='Rockfish')

  SpeciesTable<- unique(Data [,colnames(Data)[LifeStart:dim(Data)[2]]])


  SpeciesSightings<- Data %>%
    group_by(Site,CommName) %>%
    summarize(There=length(CommName))

  SpeciesSightingsByTrip<- Data %>%
    group_by(sample_Idcellday) %>%
    summarize(SpeciesSeen=length(unique(CommName)))


  BlankVars<- colnames(Data)

  BlankVars<- BlankVars[!(BlankVars%in%c('Year','Month','Site','sample_Idcellday','Sample_Type','Sample_Area','Area_units',
                                         'Angler_hours','MPA_or_REF','GRID_ID_cell','Meters.to.MPA.border','MeanLon','MeanLat','SiteId'))]
  Sites<- unique(Data$Site)

  for (s in 1:length(Sites))
  {

    Trips<- unique(Data$sample_Idcellday[Data$Site==Sites[s]])

    SpeciesSeen<- SpeciesSightings$CommName[SpeciesSightings$Site==Sites[s]]

    for (t in 1:length(Trips))
    {

      SpeciesSpotted<- unique(Data$CommName[Data$sample_Idcellday==Trips[t]])

      SpeciesMissing<- SpeciesSeen[!(SpeciesSeen %in% SpeciesSpotted)]

      if (sum(!(SpeciesSeen %in% SpeciesSpotted))>0)
      {
        BlankTrip<- Data[Data$sample_Idcellday==Trips[t],][1,]

        BlankTrip<- RepMat(BlankTrip,length(SpeciesMissing))

        BlankTrip[,colnames(BlankTrip) %in% BlankVars]<- NA

        MissingData<- SpeciesTable[SpeciesTable$CommName %in% SpeciesMissing,]

        BlankTrip[,'length_cm']<- 0

        BlankTrip [,colnames(MissingData)]<- MissingData

        Data<- rbind(Data,BlankTrip)

      } #Close if any missing loop

    } #Close Trips

  } #Close sites

  Data<- Data[order(Data$Year,Data$Month,Data$Site,Data$sample_Idcellday),]

  return(Data)
} #Close function

