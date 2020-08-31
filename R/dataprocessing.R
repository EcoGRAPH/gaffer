dataprocessing <- function(laglen   = NULL,
                           dlagdeg  = NULL,
                           nprincomps =NULL,
                           modeldata = NULL,
                           env = NULL ){

  datalagger <- expand.grid(placeid=unique(modeldata$placeid),
                            date=unique(modeldata$date),
                            lag=seq(from=0, to=laglen-1, by=1))
  #print(datalagger)
  datalagger$date-datalagger$lag
  datalagger$lag
  datalagger$date
  datalagger$laggeddate <- datalagger$date-datalagger$lag
  datalagger <- left_join(datalagger, env,
                          by=c("placeid"="placeid",
                               "laggeddate"="date"))

  envnames <- colnames(env)
  envnames <- envnames[!(envnames %in% c("placeid", "date", "year", "doy"))]
  datalagger

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

    # replace with summaries
    modeldata[,paste(curenv, "mat", sep="")] <- modeldata[,paste(curenv, "mat", sep="")] %*%
      as.matrix(bs(0:(laglen-1), intercept=TRUE))

    # and get rid of the lagged data columns
    modeldata[,grep(x=colnames(modeldata), pattern=paste(curenv, "_", sep=""), fixed=TRUE)] <- NULL

  }


  return (list(modeldata,env))

}
