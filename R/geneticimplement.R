geneticimplement <- function(individpergeneration = NULL,
                             initialclusters      = NULL,
                             initialcovars        = NULL,
                             generations          = NULL,
                             modeldata            = NULL,
                             envdata              = NULL,
                             shapefile            = NULL) {

  adjacency <- shapefile

  # retain only those shapes that appear in the model data
  adjacency <- adjacency[adjacency$placeid %in% unique(modeldata$placeid),]
  # retain a list of placeids for later
  placeids <- adjacency$placeid
  # calculate an actual adjacency list
  adjacency <- nb2listw(poly2nb(adjacency,
                                queen=TRUE,
                                row.names=adjacency$placeid),
                        style="B")

  # set up the first generation of models
  modelsdf <- data.frame(generation=1,
                         mutation="firstgen",
                         modelnumber=1:individpergeneration,
                         modelmeasure=Inf,
                         selectionprobability=1/individpergeneration)

  # set up initial clusters and covariates
  envnames <- colnames(env)
  # this next line shouldn't be hardcoded
  envnames <- envnames[!(envnames %in% c("placeid", "date", "year", "doy"))]

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

  }

  # set up initial model measures
  modelsdf$modelmeasure <- Inf
  modelsdf$selectionprobability <- 1
  modelsdf$generation <- 1

  # evaluate the first generation
  modelsdf$modelmeasure <- evaluategeneration(models=modelsdf,
                                              modeldata=mal,
                                              adjacency=adjacency,
                                              placeids=placeids)

  # rank models
  modelsdf$selectionprobability <- rank(modelsdf$modelmeasure, ties.method="random")
  modelsdf$selectionprobability <- 1/modelsdf$selectionprobability / sum(1/modelsdf$selectionprobability, na.rm=TRUE)
  modelsdf$selectionprobability[is.na(modelsdf$selectionprobability)] <- 0

  # save the first generation
  write.csv(modelsdf, ".\\outputs\\generation_1.csv")

  gens <- c()
  measure <- c()

  count = 0

  # run the rest of the generations
  for (generation in 2:generations) {

    modelsdf <- addgeneration(models=modelsdf,
                              modeldata=mal,
                              covariatenames=envnames,
                              individualspergeneration=individpergeneration,
                              adjacency=adjacency,
                              placeids=placeids)

    #if (generation == generations) {

      write.csv(modelsdf,
                paste(".\\outputs\\generation_",
                      generation, ".csv",
                      sep=""))

    #}


      gens <- c(gens, generation)
      measure <- c(measure,min(modelsdf$modelmeasure))


  }

  for (i in 1:length(gens)){

    count = count + 1
    if (count == 10){

      plot(head(gens,10),head(measure,10),type="l", main="generations vs measure",xlab = "Generations",
           ylab="Model min measure")
     # print(gens[1:2])
     # print(measure[1:2])
      gens <- gens[-seq(1:10)]
      measure <- measure[-seq(1:10)]

      count = 0
    }
  }

  return(modelsdf)

}
