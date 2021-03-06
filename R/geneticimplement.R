geneticimplement <- function(individpergeneration = NULL,
                             initialclusters      = NULL,
                             initialcovars        = NULL,
                             generations          = NULL,
                             modeldata            = NULL,
                             envnames             = NULL,
                             shapefile            = NULL,
                             restartfilename      = NULL,
                             forcecovariate       = NULL,
                             forcecyclicals       = NULL,
                             slice                = 1) {

  adjacency <- shapefile

  # retain only those shapes that appear in the model data
  adjacency <- adjacency[adjacency$placeid %in% unique(modeldata$placeid),]
  # retain a list of placeids for later
  placeids <- adjacency$placeid

  # set up covariates
  envnames <- envnames[!(envnames %in% c("placeid", "date", "year", "doy"))]

  # if we do not have a file to load, create the first generation
  if (is.null(restartfilename)) {

    # set up the first generation of models
    modelsdf <- data.frame(generation=1,
                           mutation="firstgen",
                           environdf=5,
                           longtermdf=5,
                           modelnumber=1:individpergeneration,
                           modelmeasure=Inf,
                           selectionprobability=1/individpergeneration)

    modelsdf$clusterseeds <- ""
    modelsdf$covars       <- ""
    for (i in 1:nrow(modelsdf)) {

      modelsdf$clusterseeds[i] <- paste(sample(unique(modeldata$placeid),
                                               size=initialclusters,
                                               replace=FALSE),
                                        1:initialclusters,
                                        sep=",",
                                        collapse=",")

      modelsdf$covars[i] <- paste(sample(envnames,
                                         size=initialcovars,
                                         replace=FALSE),
                                  collapse=",")

      if (is.null(forcecyclicals)) {

        modelsdf$cyclicals[i] <- sample(c("none", "percluster", "perplaceid"), size=1)

      } else {

        modelsdf$cyclicals[i] <- forcecyclicals

      }

    }
    if (!is.null(forcecovariate)) {

      modelsdf$covars <- forcecovariate

    }

    # set up initial model measures
    modelsdf$modelmeasure <- Inf
    modelsdf$selectionprobability <- 1
    modelsdf$generation <- 1

    # evaluate the first generation
    modelsdf$modelmeasure <- evaluategeneration(models=modelsdf,
                                                modeldata=mal,
                                                shapefile=shapefile,
                                                placeids=placeids)

    # rank models
    modelsdf$selectionprobability <- rank(modelsdf$modelmeasure, ties.method="random")
    modelsdf$selectionprobability <- 1/modelsdf$selectionprobability / sum(1/modelsdf$selectionprobability, na.rm=TRUE)
    modelsdf$selectionprobability[is.na(modelsdf$selectionprobability)] <- 0

    # save the first generation
    write.csv(modelsdf, ".\\csv outputs\\generation_1.csv")

    # set the minimum generation
    mingeneration <- 2

  } else {

    modelsdf <- read.csv(restartfilename)
    mingeneration <- max(modelsdf$generation, na.rm=TRUE) + 1

  }

  # run the rest of the generations
  for (generation in mingeneration:generations) {

    modelsdf <- addgeneration(models=modelsdf,
                              modeldata=mal,
                              covariatenames=envnames,
                              individualspergeneration=individpergeneration,
                              shapefile=shapefile,
                              placeids=placeids,
                              forcecovariate=forcecovariate,
                              forcecyclicals=forcecyclicals)#,
                              #slice=slice)

    if (generation %% slice == 0) {

      write.csv(modelsdf,
                paste(".\\csv outputs\\generation_",
                      generation, ".csv",
                      sep=""))

      thisplot <- ggplot(modelsdf) + geom_point(aes(x=generation,
                                                    y=rank(modelmeasure)))
      ggsave(plot=thisplot,
             filename=paste(".\\png outputs\\generation_",
                            generation, ".png",
                            sep=""))


    }

  }

  return(modelsdf)

}
