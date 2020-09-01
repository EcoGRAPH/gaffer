dataprocessing <- function(laglen = NULL,
                           modeldata = NULL,
                           env = NULL){

  datalagger <- expand.grid(placeid=unique(modeldata$placeid),
                            date=unique(modeldata$date),
                            lag=seq(from=0, to=laglen-1, by=1))

  datalagger$laggeddate <- datalagger$date-datalagger$lag
  datalagger <- left_join(datalagger, env,
                          by=c("placeid"="placeid",
                               "laggeddate"="date"))

  envnames <- colnames(env)
  envnames <- envnames[!(envnames %in% c("placeid", "date", "year", "doy"))]

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

  }

  # create a lagmat for use in regressions
  modeldata$lagmat <- matrix(rep(0:(laglen-1),
                                 nrow(modeldata)),
                             nrow(modeldata),
                             laglen,
                             byrow=TRUE)

  return(list(modeldata,env))

}
