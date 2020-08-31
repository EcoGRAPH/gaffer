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

    # get the cluster seeds
    curclusterseeds <- unlist(strsplit(x=models$clusterseeds,
                                       split=",",
                                       fixed=TRUE))
    curclusters <- data.frame(placeid=rep("", length(curclusterseeds)/2),
                              cluster=rep(0 , length(curclusterseeds)/2))
    for (i in 1:(length(curclusterseeds)/2)) {

      curclusters$placeid[i] <- curclusters[1+2*(i-1)]
      curclusters$cluster[i] <- as.numeric(curclusters[2*i])

    }
    curclusters$cluster <- factor(curclusters$cluster)

    # get the covariates
    curcovars <- unlist(strsplit(x=models$covars,
                                 split=",",
                                 fixed=TRUE))
    covarformula <- paste("s(",
                          curcovars,
                          ", by=lagmat, bs='tp')",
                          sep="",
                          collapse="+")

    # fill in the missing placeids with NAs
    curclusters <- left_join(data.frame(placeid=placeids),
                             curclusters)
    # then propagate
    curclusters$cluster <- fillbynearest(adjacency=adjacency,
                                         covariate=models$clustermat[curmodelnum,])
    modeldata <- left_join(modeldata,
                           curclusters,
                           by="placeid")

    # make sure our factors are indeed factors
    modeldata$cluster <- factor(modeldata$cluster)
    modeldata$placeid <- factor(modeldata$placeid)

    tryCatch({

      # fit the models
      modelfit <- batch_bam(data = modeldata,
                            bamargs = list("formula" = myformula1,
                                           "family" = gaussian(),
                                           "discrete" = TRUE,
                                           "nthread" = parallel::detectCores(logical=FALSE)-1),
                            bamargs_fallback = list("formula" = myfallbackformula1),
                            over = "cluster")
      #print(modelfit)
      myAICs <- extractAIC.batch_bam(models=modelfit)
      models$modelmeasure[curmodelnum] <- sum(myAICs[,2]) + (log(nrow(modeldata)) - 2)*sum(myAICs[,1])

      # try cleaning up
      rm(modelfit)
      rm(curclusters)
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
