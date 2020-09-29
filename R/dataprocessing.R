dataprocessing <- function(laglen = NULL,
                           modeldata = NULL,
                           env = NULL){

  # get list of environmental names
  envnames <- colnames(env)
  envnames <- envnames[!(envnames %in% c("placeid", "date", "year", "doy"))]

  # make sure this is a factor
  env$placeid <- factor(env$placeid)

  # complete and anomalize the data
  for (curenv in envnames) {

    env$temp <- env[,curenv]
    # thisbatch <- batch_bam(data=env,
    #                        bamargs=list("formula"=as.formula("temp ~ s(doy, bs='cc')"),
    #                                     "discrete"=TRUE,
    #                                     "nthread" = parallel::detectCores(logical=FALSE)-1),
    #                        bamargs_fallback=list("formula"=as.formula("temp ~ 1")),
    #                        over="placeid")
    #
    # env$pred <- clusterapply::predict.batch_bam(models=thisbatch,
    #                                             newdata=env,
    #                                             over="placeid")

    tempdf <- dplyr::summarise(group_by(env,
                                        placeid,
                                        doy),
                               meantemp = mean(temp, na.rm=TRUE))
    env <- left_join(env, tempdf, by=c("placeid", "doy"))
    rm(tempdf)

    # replace missing
    env[is.na(env$temp),curenv] <- env$meantemp[is.na(env$temp)]

    # anomalize
    env[,paste("anom_", curenv, sep="")] <- env[,curenv] - env$meantemp
    env[is.na(env[,paste("anom_", curenv, sep="")]),paste("anom_", curenv, sep="")] <- 0

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

  }

  # create a lagmat for use in regressions
  modeldata$lagmat <- matrix(rep(0:(laglen-1),
                                 nrow(modeldata)),
                             nrow(modeldata),
                             laglen,
                             byrow=TRUE)

  return(list(modeldata,env))

}
