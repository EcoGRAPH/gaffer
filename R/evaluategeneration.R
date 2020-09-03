evaluategeneration <- function(models=NULL,
                               modeldata=NULL,
                               adjacency=NULL,
                               placeids=NULL) {

  # make sure we start from scratch
  models$modelmeasure <- Inf
  models$selectionprobability <- 1

  # create some helpful covariates
  modeldata$numdate <- as.numeric(modeldata$date)
  modeldata$doy     <- as.numeric(format(modeldata$date, "%j"))

  # create the basic formula
  baseformula <- "objective ~ placeid + s(numdate, by=placeid, bs='tp', id=1) + s(doy, bs='cc', id=2)"
  basefallback <- "objective ~ s(numdate, id=1) + s(doy, bs='cc', id=2)"

  # run batch_bam on all the models
  for (curmodelnum in 1:nrow(models)) {

    # get rid of clusters from previous model, if any
    modeldata$cluster <- NULL

    # get the cluster seeds
    curclusterseeds <- unlist(strsplit(x=models$clusterseeds[curmodelnum],
                                       split=",",
                                       fixed=TRUE))
    curclusters <- data.frame(placeid=rep("", length(curclusterseeds)/2),
                              cluster=rep(0 , length(curclusterseeds)/2))
    for (i in 1:(length(curclusterseeds)/2)) {

      curclusters$placeid[i] <- curclusterseeds[2*i-1]
      curclusters$cluster[i] <- as.numeric(curclusterseeds[2*i])

    }
    curclusters$cluster <- factor(curclusters$cluster)

    # get the covariates
    curcovars <- unlist(strsplit(x=models$covars[curmodelnum],
                                 split=",",
                                 fixed=TRUE))
    if (curcovars != "none") {

      covarformula <- paste("s(",
                            curcovars,
                            "mat, by=lagmat, bs='tp')",
                            sep="",
                            collapse="+")
      modelformula <- as.formula(paste(baseformula, covarformula, sep="+"))

    } else {

      modelformula <- as.formula(baseformula)

    }

    # create the formulas
    fallbackformula <- as.formula(basefallback)

    # fill in the missing placeids with NAs
    curclusters <- left_join(data.frame(placeid=placeids),
                             curclusters,
                             by="placeid")

    # then propagate
    curclusters$cluster <- fillbynearest(adjacency=adjacency,
                                         covariate=curclusters$cluster)
    modeldata <- left_join(modeldata,
                           curclusters,
                           by="placeid")

    # make sure our factors are indeed factors
    modeldata$cluster <- factor(modeldata$cluster)
    modeldata$placeid <- factor(modeldata$placeid)

    tryCatch({

      # fit the models
      modelfit <- batch_bam(data = modeldata,
                            bamargs = list("formula" = modelformula,
                                           "family" = gaussian(),
                                           "discrete" = TRUE,
                                           "nthread" = parallel::detectCores(logical=FALSE)-1),
                            bamargs_fallback = list("formula" = fallbackformula),
                            over = "cluster")
      myAICs <- extractAIC.batch_bam(models=modelfit)
      models$modelmeasure[curmodelnum] <- sum(myAICs[,2]) + (2*log(nrow(modeldata)) - 2)*sum(myAICs[,1])
      #models$modelmeasure[curmodelnum] <- sum(myAICs[,2])

      #save(modelfit, file="modelfit.rdata")

      # try cleaning up
      rm(modelfit)
      rm(curclusters)
      rm(modelformula)
      rm(fallbackformula)
      modeldata$cluster <- NULL
      gc()

    }, error = function(e) {

      models$modelmeasure[curmodelnum] <- Inf
      modeldata$cluster <- NULL
      rm(list=ls())
      gc()
      return(NULL)

    })

  }

  # retain smallest possible frame
  returnframe <- unlist(models$modelmeasure)
  rm(list=setdiff(ls(), "returnframe"))
  gc()

  return(returnframe)

}
