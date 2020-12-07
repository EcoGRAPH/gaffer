dataprocessing <- function(laglen = NULL,
                           modeldata = NULL,
                           env = NULL,
                           maxdeg = 10){

  # get list of environmental names
  envnames <- colnames(env)
  envnames <- envnames[!(envnames %in% c("placeid", "date", "year", "doy"))]

  # make sure this is a factor
  env$placeid <- factor(env$placeid)

  # complete and anomalize the data
  for (curenv in envnames) {

    env$temp <- env[,curenv]
    tempdf <- dplyr::summarise(group_by(env,
                                        placeid,
                                        doy),
                               meantemp = mean(temp, na.rm=TRUE))
    env <- left_join(env, tempdf, by=c("placeid", "doy"))
    rm(tempdf)

    ## replace missing
    env[is.na(env$temp),curenv] <- env$meantemp[is.na(env$temp)]
    env[is.na(env$temp),curenv] <- mean(env$temp, na.rm=TRUE)

    # anomalize
    env[,curenv] <- env[,curenv] - env$meantemp
    env[is.na(env[,curenv]),curenv] <- 0

    # remove temporary variables
    env$meantemp <- NULL
    env$temp     <- NULL

  }

  # get list of environmental names again, now with anoms
  envnames <- colnames(env)
  envnames <- envnames[!(envnames %in% c("placeid", "date", "year", "doy"))]

  datalagger <- expand.grid(placeid=unique(modeldata$placeid),
                            date=unique(modeldata$date),
                            lag=seq(from=0, to=laglen-1, by=1))

  datalagger$laggeddate <- datalagger$date-datalagger$lag
  datalagger <- left_join(datalagger, env,
                          by=c("placeid"="placeid",
                               "laggeddate"="date"))

  for (curenv in envnames) {

    # pivot environmental data by lag
    envlagged <- dcast(datalagger, placeid + date ~ lag, value.var=curenv)
    names(envlagged) <- paste(curenv,names(envlagged),sep="_")

    # join to malaria data
    modeldata <- left_join(x=modeldata, y=envlagged,
                           by=c("placeid"=paste(curenv,"placeid",sep="_"),
                                "date"=paste(curenv,"date",sep="_")))

    # turn these lagged data into matrices
    modeldata$tempvar <- as.matrix(modeldata[,grep(x=colnames(modeldata),
                                                   pattern=paste(curenv, "_", sep=""),
                                                   fixed=TRUE)])
    colnames(modeldata)[length(colnames(modeldata))] <- paste(curenv, "mat", sep="")

    # and get rid of the lagged data columns
    modeldata[,grep(x=colnames(modeldata), pattern=paste(curenv, "_", sep=""), fixed=TRUE)] <- NULL

    # create the lagged summary statistics
    polybas <- poly(0:(laglen-1),
                    degree=maxdeg,
                    simple=TRUE)
    polybas <- cbind(rep(1, dim(polybas)[2]), polybas)
    modeldata[,paste(curenv, "mat", sep="")] <- modeldata[,paste(curenv, "mat", sep="")] %*% polybas

  }

  # create the cyclical summaries
  modeldata$doy <- as.numeric(format(modeldata$date, "%j"))
  for (curtotdf in 4:maxdeg) {

    modeldata[,paste("cyclicalmat_", curtotdf, sep="")] <- cSplineDes(x=modeldata$doy,
                                                                      knots=seq(from=0,
                                                                                to=366,
                                                                                length.out=curtotdf+1))
  }

  # # create a lagmat for use in regressions
  # modeldata$lagmat <- matrix(rep(0:(laglen-1),
  #                                nrow(modeldata)),
  #                            nrow(modeldata),
  #                            laglen,
  #                            byrow=TRUE)

  return(list(modeldata,env))

}
