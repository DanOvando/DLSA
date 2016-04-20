CatchMSY<- function(Data,n,ProcessError,phi = 0.188,resilience = 'high',Smooth,
                    Display = F,CatchJitters,CatchError,CatchBias,InterpCatch,
                    StartYear,startbio = c(0.6,1),MidYear,midbio = c(0,1),finalbio)
{
  ### Run CatchMSY, adapted form Martell & Froese 2012
  # CatchDat=CatchData
  # n=5000
  # ProcessError=0.05
  # Smooth=0
  # Display=1
  # CatchJitters=1
  # CatchError=0.075
  # CatchBias= 0 #the mean of the log normal errors in catch, 0 means no bias, above 0 positive bias, - negative bias
  # InterpCatch=1
  # StartYear=1994
  # StartBio=c(0.5,0.75)
  # MidYear=NA
  # MidBio=NA
  # FinalBio=c(0.25,0.5)
  # # Fish$res='Low'
  # CatchDat: Time series of catch
  # n: The number of iterations to search over
  # ErrorSize: The maginitude of the gradient of life history terms to search over
  # Smooth: Marks whether or not to smooth catch history
  # Display: Show diagnostics as it goes
  # CatchJitters: Number of times to try modifications of the catch
  # CatchError: Log error magnitude


  Data<- if (is.na(StartYear)){Data} else {Data[Data$Year>=StartYear,]}

  yr   <- Data$Year

  set.seed(999)  ## for same random sequence

  Output<- as.data.frame(matrix(NA,nrow=length(yr),ncol=9))

  colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')

  MCOutput<- as.data.frame(matrix(NA,nrow=length(yr),ncol=7))

  colnames(MCOutput)<- c('Iteration','Year','Method','SampleSize','Value','Metric','Flag')

  MCDetails<- as.data.frame(matrix(NA,nrow=length(yr)*CatchJitters,ncol=2))

  colnames(MCDetails)<- c('Year','FvFmsy')

  Data$RanCatchMSY<- FALSE

  Data$phi <- phi

  Data$resilience <- resilience

  Data$BtoKRatio<- 1/((Data$phi+1)^(1/Data$phi))

  Data$g<- NA

  Data$k<- NA

  Data$MSYLogSd<- NA

  Data$gLogSd<- NA

  Data$KLogSd<- NA

  Data$CatchMSYBvBmsy<- NA

  Data$CatchMSYBvBmsy_LogSd<- NA

  # Run Catch MSY -----------------------------------------------------------

  require(zoo, quietly)

  MatrixCmsy<- function(parbound,n,interbio,finalbio,startbt)
  {

    with(as.list(parbound),
         {

           gi = rep(exp(runif(n, log(start_g[1]), log(start_g[2]))),length(startbt))  ## get N values between g[1] and g[2], assign to ri

           ki = rep(exp(runif(n, log(start_k[1]), log(start_k[2]))),length(startbt))  ## get N

           startbti<- sort(rep(startbt,n))

           ParamSpace<- as.data.frame(cbind(phi,gi,ki,interbio[1],interbio[2],finalbio[1],finalbio[2],sigR,startbti))

           colnames(ParamSpace)<- c('phi','g','K','InterBio1','InterBio2','FinalBio1','FinalBio2','sigR','StartBio')

           #     bvbmsyell<- rep(0,length(bt))
           #     if(bt[nyr]/k>=lam1 && bt[nyr]/k <=lam2 && min(bt,na.rm=T) > 0 && max(bt,na.rm=T) <=k && bt[which(yr==interyr)]/k>=interbio[1] && bt[which(yr==interyr)]/k<=interbio[2])
           #     {ell = 1
           #      bvbmsyell<-rep(1,length(bt))
           #     }

           CatchMat<- matrix(rep(ct,dim(ParamSpace)[1]),nrow=dim(ParamSpace)[1],ncol=length(ct),byrow=T)

           btMat<- matrix(NA,nrow=dim(CatchMat)[1],dim(CatchMat)[2])

           btMat[,1]<- ParamSpace$K*ParamSpace$StartBio*exp(rnorm(dim(btMat)[1],0, ParamSpace$sigR))

           for (y in 2:length(ct))
           {
             xt <- exp(rnorm(dim(btMat)[1],0, sigR))

             btMat[,y]<- (btMat[,y-1]+((phi+1)/phi)*ParamSpace$g*btMat[,y-1]*(1-(btMat[,y-1]/ParamSpace$K)^phi)-ct[y-1])*(xt)
           }

           ItId<- 1:dim(btMat)[1]

           btMat <- as.data.frame(btMat)

           colnames(btMat) <- paste('X',yr,sep = '')

           ResultMat<- data.frame(ItId,btMat,ParamSpace)

           BioDat<- ResultMat[,grepl('X',colnames(ResultMat))]

           interyr<- round(median(1:nyr))

           EllBio<- data.frame(apply(BioDat,1,min),apply(BioDat,1,max),BioDat[,interyr]/ResultMat$K,BioDat[,nyr]/ResultMat$K)

           colnames(EllBio)<- c('MinBio','MaxBio','InterBio','FinalBio')

           #   Ell= EllBio$FinalBio>=ResultMat$FinalBio1 & EllBio$FinalBio <= ResultMat$FinalBio2 & EllBio$InterBio>=ResultMat$InterBio1 & EllBio$InterBio <= ResultMat$InterBio2 & EllBio$MinBio>0 & EllBio$MaxBio<ResultMat$K

           Ell= ResultMat$StartBio==min(ResultMat$StartBio) & EllBio$FinalBio>=ResultMat$FinalBio1 & EllBio$FinalBio <= ResultMat$FinalBio2 & EllBio$InterBio>=ResultMat$InterBio1 & EllBio$InterBio <= ResultMat$InterBio2 & EllBio$MinBio>0 & EllBio$MaxBio<ResultMat$K

           Missing<- is.na(EllBio$FinalBio)

           PossibleRuns<- ResultMat[Ell & !Missing,]
           return(PossibleRuns)

         })
  }

  BtoKRatio<- unique(Data$BtoKRatio)

  CatchYears<- (Data$Year*as.numeric(is.na(Data$Catch)==F))

  CatchYears[CatchYears==0]<- NA

  FirstCatchYear<- which(Data$Year==min(CatchYears,na.rm=T))[1]

  LastCatchYear<- which(Data$Year==max(CatchYears,na.rm=T))[1]

  Data<- Data[FirstCatchYear:LastCatchYear,]

  yr   <- Data$Year

  ct   <- Data$Catch  ## assumes that catch is given in tonnes, transforms to 1'000 tonnes

  PossibleRuns<- NA

  if (sum(ct,na.rm=T)>0 & length(LastCatchYear)>0 & length(ct)>1)
  {

    ct<- na.approx(ct)

    if(Smooth==1){ct<- runmed(ct,3)}

    res  <- unique(Data$resilience)

    if(is.na(res)){res<- 0.5}

    for (i in 1){
      start_g  <- if(res == "Very low"){c(0.001, 0.05)}
      else if(res == "Low") {c(0.05,0.15)}
      else if(res == "Medium") {c(0.15,0.5)}
      else if(res == "High") {c(0.5,1)}
      else {c(0.15,0.5)} ## Medium, or default if no res is found
    }

    phi<- unique(Data$phi)

    start_g<- start_g*(phi/(1+phi)) #To account for g instead of r

    nyr  <- length(yr)    ## number of years in the time series

    flush.console()

    ## PARAMETER SECTION

    start_k     <- c(max(ct,na.rm=T),50*max(ct,na.rm=T)) ## default for upper k e.g. 100 * max catch

    if (all(is.na(startbio)))
    {
      startbio    <- if(ct[1]/max(ct,na.rm=T) < 0.5) {c(0.5,0.9)} else {c(0.3,0.6)} ## use for batch processing #SUB IN BVBMSY VALUES
    }

    interyr 	<- median(1:length(yr))   ## interim year within time series for which biomass estimate is available; set to yr[2] if no estimates are available #SUB IN INTERMIN YEAR

    interbio   <- midbio ## biomass range for interim year, as fraction of k; set to 0 and 1 if not available

    interyr<- yr[interyr]

    if (all(is.na(finalbio)))
    {
      warning('Using automated final depletion: you get what you put in!')
      finalbio    <- if(ct[nyr]/max(ct,na.rm=T) > 0.5) {c(0.3,0.7)} else {c(0.01,0.4)} ## use for batch processing #SET TO KNOWN B/BMSY RANGE
    }

    sigR <- ProcessError
    #       sigR        <- 0.0      ## process error; 0 if deterministic model; 0.05 reasonable value? 0.2 is too high
    #     startbt     <- seq(startbio[1], startbio[2], length.out = 10) ## apply range of start biomass in steps of 0.05
    startbt     <- seq(startbio[1], startbio[2], by = 0.05) ## apply range of start biomass in steps of 0.05

    parbound <- list(g = start_g, k = start_k, lambda = finalbio, sigR=sigR,phi=unique(Data$phi))

    if (Display == T)
    {
      cat("Last year =",max(yr),", last catch =",ct[nyr],"\n")
      cat("Resilience =",res,"\n")
      cat("Process error =", sigR,"\n")
      cat("Assumed initial biomass (B/k) =", startbio[1],"-", startbio[2], " k","\n")
      cat("Assumed intermediate biomass (B/k) in", interyr, " =", interbio[1],"-",interbio[2]," k","\n")
      cat("Assumed final biomass (B/k) =", parbound$lambda[1],"-",parbound$lambda[2]," k","\n")
      cat("Initial bounds for g =", parbound$g[1], "-", parbound$g[2],"\n")
      cat("Initial bounds for k =", format(parbound$k[1], digits=3), "-", format(parbound$k[2],digits=3),"\n")
    }
    flush.console()

    ## MAIN

    PossibleRuns <- MatrixCmsy(parbound,n,interbio,finalbio,startbt)

    ## Get statistics on g, k, MSY and determine new bounds for g and k
    g1 	<- PossibleRuns$g
    k1 	<- PossibleRuns$K

    if(length(g1)<10)
    {

      finalbio<- pmax(0,pmin(1,finalbio+c(-.065,.065)))
      PossibleRuns<- MatrixCmsy(parbound,n,interbio,finalbio,startbt)
      ## Get statistics on g, k, MSY and determine new bounds for g and k
      g1   <- PossibleRuns$g
      k1 	<- PossibleRuns$K

    }

    if(length(g1)<10) {
      cat("Too few (", length(g1), ") possible g-k combinations, check input parameters","\n")
      flush.console()
    }

    if(length(g1)>=10) {

      msy1  <- (g1*k1)*BtoKRatio
      mean_msy1 <- exp(mean(log(msy1)))
      max_k1a  <- min(k1[g1<1.1*parbound$g[1]],na.rm=T) ## smallest k1 near initial lower bound of g
      max_k1b  <- max(k1[(g1*k1)*BtoKRatio <mean_msy1],na.rm=T) ## largest k1 that gives mean MSY
      max_k1 <- if(max_k1a < max_k1b) {max_k1a} else {max_k1b}
      ## set new upper bound of g to 1.2 max r1
      parbound$g[2] <- 1.2*max(g1)
      ## set new lower bound for k to 0.9 min k1 and upper bound to max_k1
      parbound$k 	  <- c(0.9 * min(k1), max_k1)

      if (Display == T)
      {
        cat("First MSY =", format(mean_msy1, digits=3),"\n")
        cat("First g =", format(exp(mean(log(g1))), digits=3),"\n")
        cat("New upper bound for g =", format(parbound$g[2],digits=2),"\n")
        cat("New range for k =", format(parbound$k[1], digits=3), "-", format(parbound$k[2],digits=3),"\n")
      }

      ## Repeat analysis with new g-k bounds
      PossibleRuns<- MatrixCmsy(parbound,n,interbio,finalbio,startbt)

      PossibleRuns$Fail<- 0

      ## Get statistics on g, k and msy
      g   <- PossibleRuns$g
      k 	<- PossibleRuns$K

      PossibleRuns$MSY<- (g*k)*BtoKRatio


      PossibleRuns <- PossibleRuns %>%
        gather('Year','biomass', which(grepl('X',colnames(PossibleRuns)))) %>%
        mutate(Year = sub('X','',Year),bvbmsy = (biomass/k)/BtoKRatio) %>%
        mutate(Year = as.numeric(Year)) %>%
        left_join(Data[,c('Year','Catch')], by = 'Year') %>%
        mutate(fvfmsy = (Catch/MSY)/bvbmsy)

      bvbmsy <- PossibleRuns %>%
        group_by(Year) %>%
        summarize(meanb = exp(mean(log(bvbmsy))), sdb = (sd(log(bvbmsy))) )

      #       bvbmsy<- (PossibleRuns[,grepl('X',colnames(PossibleRuns))]/k)/BtoKRatio


      time_bvbmsy<- bvbmsy$meanb
      #       mean_bvbmsy<- mean(apply(bvbmsy,1,function(x) exp(mean(log(x)))))
      LogSD_bvbmsy<- bvbmsy$sdb

      msy = (g * k) * BtoKRatio

      Fmsy<- g

      mean_ln_msy = mean(log(msy),na.rm=T)

      mean_ln_g<- mean(log(g),na.rm=T)

      mean_ln_k<- mean(log(k),na.rm=T)

      Data$RanCatchMSY<- TRUE

      Data$MSY<- exp(mean_ln_msy)

      Data$g<- exp(mean_ln_g)

      Data$k <- exp(mean_ln_k)

      Data$MSYLogSd <- (sd(log(msy)))

      Data$gLogSd <- (sd(log(g),na.rm=T))

      Data$KLogSd <- (sd(log(k),na.rm=T))

      Data$CatchMSYBvBmsy <- time_bvbmsy

      Data$BvBmsy <- time_bvbmsy

      Data$CatchMSYBvBmsy_LogSd <- LogSD_bvbmsy

      Data$FvFmsy <- (Data$Catch / Data$MSY) / Data$BvBmsy

      ## plot MSY over catch data
      if (Display == T  & length(g)>10)
      {
        pdf(file = paste(FigureFolder,'CMSY Diagnostics.pdf', sep = ''))
        par(mfcol=c(2,3))
        plot(yr, ct, type="l", ylim = c(0, max(ct)), xlab = "Year", ylab = "Catch (MT)", main = unique(Data$Species))
        abline(h=exp(mean(log(msy))),col="red", lwd=2)
        abline(h=exp(mean_ln_msy - 2 * sd(log(msy))),col="red")
        abline(h=exp(mean_ln_msy + 2 * sd(log(msy))),col="red")
        hist(g, freq=F, xlim=c(0, 1.2 * max(g,na.rm=T)), main = "")
        abline(v=exp(mean(log(g))),col="red",lwd=2)
        abline(v=exp(mean(log(g))-2*sd(log(g))),col="red")
        abline(v=exp(mean(log(g))+2*sd(log(g))),col="red")

        plot(g1, k1, xlim = start_g, ylim = start_k, xlab="g", ylab="k (MT)")
        hist(k, freq=F, xlim=c(0, 1.2 * max(k)), xlab="k (MT)", main = "")
        abline(v=exp(mean(log(k))),col="red", lwd=2)
        abline(v=exp(mean(log(k))-2*sd(log(k))),col="red")
        abline(v=exp(mean(log(k))+2*sd(log(k))),col="red")

        plot(log(g), log(k),xlab="ln(g)",ylab="ln(k)")
        abline(v=mean(log(g)))
        abline(h=mean(log(k)))
        abline(mean(log(msy))+log(4),-1, col="red",lwd=2)
        abline(mean(log(msy))-2*sd(log(msy))+log(4),-1, col="red")
        abline(mean(log(msy))+2*sd(log(msy))+log(4),-1, col="red")

        hist(msy, freq=F, xlim=c(0, 1.2 * max(c(msy))), xlab="MSY (MT)",main = "")
        abline(v=exp(mean(log(msy))),col="red", lwd=2)
        abline(v=exp(mean_ln_msy - 2 * sd(log(msy))),col="red")
        abline(v=exp(mean_ln_msy + 2 * sd(log(msy))),col="red")
      }
      if (Display==T)
      {
        cat("Possible combinations = ", length(g),"\n")
        cat("geom. mean g =", format(exp(mean(log(g))),digits=3), "\n")
        cat("g +/- 2 SD =", format(exp(mean(log(g))-2*sd(log(g))),digits=3),"-",format(exp(mean(log(g))+2*sd(log(g))),digits=3), "\n")
        cat("geom. mean k =", format(exp(mean(log(k))),digits=3), "\n")
        cat("k +/- 2 SD =", format(exp(mean(log(k))-2*sd(log(k))),digits=3),"-",format(exp(mean(log(k))+2*sd(log(k))),digits=3), "\n")
        cat("geom. mean MSY =", format(exp(mean(log(msy))),digits=3),"\n")
        cat("MSY +/- 2 SD =", format(exp(mean_ln_msy - 2 * sd(log(msy))),digits=3), "-", format(exp(mean_ln_msy + 2 * sd(log(msy))),digits=3), "\n")

      }
      dev.off()
      RanCMSY<- TRUE

    } #Close if r1 is greater than 10

  } #Close if there is catch loop


  CMSYResults <- list(CatchMSY = Data, PossibleParams = PossibleRuns)

  CmsyStore<- as.data.frame(matrix(NA, nrow = 0, ncol = dim(CMSYResults$CatchMSY)[2]))

  PossibleParams <- CMSYResults$PossibleParams
  #   EmptyParams <- lapply(seq(along = PossibleParams), function(i)   sum(is.na(PossibleParams[[i]]))==0)

  #   HasData<- ldply(EmptyParams)

  #   PossibleParams<- PossibleParams[which(HasData==T)]

  CmsyStore <- CMSYResults$CatchMSY
  #   PossibleParams<- ldply(PossibleParams)
  #   if (dim(PossibleParams)[1]>0 & sum(PossibleParams$Fail==0,na.rm=T)>=1)
  #   {
  #     PossibleParams<- PossibleParams[,c('Yg','phi','K','MSY','fvfmsy','bvbmsy')]
  #   }

  ConCatDat<- paste(Data$Year,sep='-')

  ConCatCmsy<- paste(CmsyStore$Year,sep='-')

  Where<- ConCatDat %in% ConCatCmsy

  Data[Where,]<- CmsyStore

  # Process MSY Results -----------------------------------------------------

  Output<- as.data.frame(matrix(NA,nrow=length(yr),ncol=9))

  colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')

  Output$Year<- yr
  Output$Method<-'CatchMSY'
  Output$SampleSize<-ct
  Output$Metric<- 'Catch/MSY'
  Output$Value<- ct/exp(mean_ln_msy)
  Output$UpperCI<- ct/quantile(msy,0.05)
  Output$LowerCI<- ct/quantile(msy,0.95)
  Output$SD<- sd((msy))

  catch_over_msy_plot <- (ggplot(PossibleParams, aes(x = factor(Year), y = Catch/MSY)) +
                            geom_violin(aes(fill = Catch)) +
                            scale_fill_continuous(low = 'green', high = 'red') +
                            geom_hline(aes(yintercept = 1)) +
                            xlab('Year') +
                            ylab('Captura / RMS') +
                            Theme + theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0.9, vjust = 0.9 )))
  ggsave(file = paste(FigureFolder, 'CatchMSY catch over MSY.pdf', sep = ''), plot = catch_over_msy_plot)

  time_trend <- PossibleParams %>%
    group_by(Year) %>%
    summarize(meanb = mean(bvbmsy,na.rm=T), meanf = mean(fvfmsy, na.rm = T), meancatch = mean(Catch, na.rm=T))

    cmsy_kobe_plot <- (ggplot(time_trend, aes(x = meanb, y = meanf)) +
                       geom_smooth(se = F, size = 1.2, color = 'black', alpha = 0.2) +
                       geom_point(aes(fill = Year, size = meancatch), pch = 21) +
                       scale_fill_gradient(low = 'blue', high = 'green')+
                       geom_text(aes(label = Year), hjust = -0.1, alpha = 0.5) +
                       geom_hline(aes(yintercept = 1)) +
                       geom_vline(aes(xintercept = 1)))
  ggsave(file = paste(FigureFolder,'CMSY Kobe Plot.pdf', sep = ''), plot = cmsy_kobe_plot, height = 6, width = 8)




  #   MCDetails$Year<- yr
  #
  #   MCDetails$FvFmsy<- Data$FvFmsy
  #
  #   MCDetails$LowerFvFmsy<- FvFmsyBox$stats[1,]
  #
  #   MCDetails$UpperFvFmsy<- FvFmsyBox$stats[5,]
  #
  #
  #   MCDetails$BvBmsy<- BvBmsyBox$stats[3,1:nyr]
  #
  #   MCDetails$LowerBvBmsy<- BvBmsyBox$stats[1,1:nyr]
  #
  #   MCDetails$UpperBvBmsy<- BvBmsyBox$stats[5,1:nyr]
  #
  #   ## Write results into outfile, in append mode (no header in file, existing files will be continued)
  #   output = data.frame( sigR, startbio[1], startbio[2], interbio[1], interbio[2], finalbio[1], finalbio[2], min(yr), max(yr), res, max(ct), ct[1], ct[nyr], length(r), exp(mean(log(r))), sd(log(r)), min(r), quantile(r,0.05), quantile(r,0.25), median(r), quantile(r,0.75), quantile(r,0.95), max(r), exp(mean(log(k))), sd(log(k)), min(k), quantile(k, 0.05), quantile(k, 0.25), median(k), quantile(k, 0.75), quantile(k, 0.95), max(k), exp(mean(log(msy))), sd(log(msy)), min(msy), quantile(msy, 0.05), quantile(msy, 0.25), median(msy), quantile(msy, 0.75), quantile(msy, 0.95), max(msy))
  #
  #   write.csv(output, file = paste(ResultFolder,'Raw CatchMSY Output.csv',sep=''), row.names = FALSE)
  #


  return(Results = list(Output = Output, PossibleParams = PossibleParams))

}




