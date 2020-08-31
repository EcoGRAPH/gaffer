geneticimplement <- function(individpergeneration = NULL,
                             initialclusters      = NULL,
                             initialcovars        = NULL,
                             generations          = NULL,
                             modeldata            = NULL,
                             envdata              = NULL,
                             shapefile            = NULL,
                             numofcovnts          = NULL) {

  adjacency <- shapefile

  # retain only those shapes that appear in the model data
  # adjacencey newpcode shouldn't be hardcoded - needs fixes throughout
  adjacency <- adjacency[adjacency$NewPCODE %in% unique(modeldata$placeid),]
  adjacency <- nb2listw(poly2nb(adjacency,
                                queen=TRUE,
                                row.names=adjacency$NewPCODE),
                        style="B")

  # set up the first generation of models
  modelsdf <- data.frame(generation=1,
                         mutation="firstgen",
                         modelnumber=1:individpergeneration,
                         modelmeasure=-Inf,
                         selectionprobability=1/individpergeneration)

  # set up initial clusters and covariates
  envnames <- colnames(env)
  # this next line shouldn't be hardcoded
  envnames <- envnames[!(envnames %in% c("placeid", "date", "year", "doy"))]

  modelsdf$clusterseeds <- ""
  modelsdf$covars       <- ""
  for (i in 1:nrow(modelsdf)) {

    modelsdf$clusterseeds[i] <- paste(sample(modeldata$placeid,
                                             size=initialclusters,
                                             replace=FALSE),
                                      1:initialclusters,
                                      sep=",",
                                      collapse=",")

    modelsdf$covars[i] <- paste(sample(envnames,
                                       size=initialcovars,
                                       replace=FALSE),
                                collapse=",")

  }

  modelsdf$modelmeasure <- Inf
  modelsdf$selectionprobability <- 1
  modelsdf$generation <- 1

  # evaluate the first generation
  modelsdf$modelmeasure <- evaluategeneration(models=modelsdf,
                                              modeldata=mal,
                                              adjacency=adjacency)

  # rank models
  modelsdf$selectionprobability <- rank(modelsdf$modelmeasure, ties.method="random")
  modelsdf$selectionprobability <- 1/modelsdf$selectionprobability / sum(1/modelsdf$selectionprobability, na.rm=TRUE)
  modelsdf$selectionprobability[is.na(modelsdf$selectionprobability)] <- 0

  # run the rest of the generations
  for (generation in 2:generations) {

    modelsdf <- addgeneration(models=modelsdf,
                               modeldata=mal,
                               covariatenames=envnames,
                               individualspergeneration=individpergeneration,
                               adjacency=adjacency,
                               numofcovnts = numofcovnts)

    if (generation == generations) {

      write.csv(modelsdf %>% dplyr::select(!starts_with("clustermat", ignore.case = TRUE)),
                paste("generation_",
                      generation, ".csv",
                      sep=""))

    }

  }

}
