evaluategeneration <- function(models=NULL,
                               modeldata=NULL,
                               shapefile=NULL,
                               placeids=NULL) {

  # make sure we start from scratch
  models$modelmeasure <- Inf
  models$selectionprobability <- 1

  # create some helpful covariates
  modeldata$numdate <- as.numeric(modeldata$date)
  modeldata$doy     <- as.numeric(format(modeldata$date, "%j"))

  # create the basic formula
  basefallback <- "objective ~ poly(numdate, degree=4)"

  # run batch_bam on all the models
  for (curmodelnum in 1:nrow(models)) {

    if (models$cyclicals[curmodelnum] == "none") {

      baseformula <- "objective ~ placeid + poly(numdate, degree=4)*placeid"

    }
    if (models$cyclicals[curmodelnum] == "percluster") {

      baseformula <- "objective ~ placeid + poly(numdate, degree=4)*placeid + cyclicalmat_5"

    }
    if (models$cyclicals[curmodelnum] == "perplaceid") {

      baseformula <- "objective ~ placeid + poly(numdate, degree=4)*placeid + cyclicalmat_5*placeid"

    }

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

      # covarformula <- paste("s(lagmat, by=",
      #                       curcovars,
      #                       "mat, bs='tp')",
      #                       sep="",
      #                       collapse="+")
      covarformula <- paste(curcovars, "mat", sep="", collapse="+")
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
    curclusters$cluster <- fillbynearest(shapefile=shapefile,
                                         covariate=curclusters$cluster)
    modeldata <- left_join(modeldata,
                           curclusters,
                           by="placeid")

    # make sure our factors are indeed factors
    modeldata$cluster <- factor(modeldata$cluster)
    modeldata$placeid <- factor(modeldata$placeid)

    # have a row number in there for stratification
    modeldata$reserved_rownum <- 1:nrow(modeldata)

    tryCatch({

<<<<<<< HEAD
      modelmeasures <- c()
      for (curtrial in 1:3) {

        # select test and training
        trainingselector <- splitstackshape::stratified(modeldata[c("placeid","reserved_rownum")],
                                                        group="placeid",
                                                        size=0.25)
        training <- modeldata[modeldata$reserved_rownum %in% trainingselector$reserved_rownum,]
        test <- modeldata[!(modeldata$reserved_rownum %in% trainingselector$reserved_rownum),]

        modelfit <- batch_lm(data = training,
                             lmargs = list("formula" = modelformula),
                             over = "cluster")
        test$preds <- predict.batch_lm(models=modelfit,
                                       over="cluster",
                                       newdata=test)

        modelmeasures[length(modelmeasures)+1] <- mean(abs(test$preds - test$objective), na.rm=TRUE)

      }

      models$modelmeasure[curmodelnum] <- mean(modelmeasures)

      # # fit the models
      # modelfit <- batch_bam(data = modeldata,
      #                       bamargs = list("formula" = modelformula,
      #                                      "family" = gaussian(),
      #                                      "discrete" = TRUE,
      #                                      "nthread" = parallel::detectCores(logical=FALSE)-1),
      #                       bamargs_fallback = list("formula" = fallbackformula),
      #                       over = "cluster")
      # myAICs <- extractAIC.batch_bam(models=modelfit)
      #
      # numclust <- max(as.numeric(modeldata$cluster), na.rm=TRUE)
      # models$modelmeasure[curmodelnum] <- sum(myAICs[,2]) + (log(nrow(modeldata)) - 2)*sum(myAICs[,1]) + 200*numclust
      # #models$modelmeasure[curmodelnum] <- sum(myAICs[,2])
=======
      # fit the models
      modelfit <- batch_bam(data = modeldata,
                            bamargs = list("formula" = modelformula,
                                           "family" = gaussian(),
                                           "discrete" = TRUE,
                                           "nthread" = parallel::detectCores(logical=FALSE)-1),
                            bamargs_fallback = list("formula" = fallbackformula),
                            over = "cluster")
      myAICs <- extractAIC.batch_bam(models=modelfit)
      #models$modelmeasure[curmodelnum] <- sum(myAICs[,2]) + (2*log(nrow(modeldata)) - 2)*sum(myAICs[,1])
      models$modelmeasure[curmodelnum] <- sum(myAICs[,2])
>>>>>>> parent of 8f94e0a... commit before RRSV

      # try cleaning up
      if (curmodelnum < nrow(models)) {

        rm(modelfit)

      }
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

  # # save the last model
  # if (savebest) {
  #
  #   save(modelfit, file=".\\model outputs\\modelfit.rdata")
  #
  # }

  # retain smallest possible frame
  returnframe <- unlist(models$modelmeasure)
  rm(list=setdiff(ls(), "returnframe"))
  gc()

  return(returnframe)

}
