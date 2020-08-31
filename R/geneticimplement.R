geneticimplement <- function(individpergeneration = NULL,
                             initialclusters      = NULL,
                             generations          =  NULL,
                             modeldata = NULL,
                             envdata = NULL,
                             shapefile = NULL,
                             whichregion = NULL,
                             numofcovnts = NULL ) {

  adjacency <- shapefile

  if (!is.null(whichregion)) {

    adjacency <- adjacency[adjacency$R_NAME == whichregion,]

  }

  adjacency2 <- adjacency
  adjacency <- adjacency[adjacency$NewPCODE %in% modeldata$placeid,]
  savenames <- adjacency$NewPCODE
  adjacency <- nb2listw(poly2nb(adjacency, queen=TRUE, row.names=adjacency$NewPCODE),
                        style="B")

  # set up the first generation of models
  modelsdf1 <- data.frame(generation=1,
                          mutation="firstgen",
                          modelnumber=1:individpergeneration,
                          modelmeasure=-Inf,
                          selectionprobability=1/individpergeneration)
  # set up initial clusters
  modelsdf1$clustermat <- matrix(rep(NA,
                                     size=nrow(modelsdf1)*length(unique(mal$placeid)),
                                     replace=TRUE),
                                 nrow=nrow(modelsdf1),
                                 ncol=length(unique(mal$placeid)))

  for (i in 1:nrow(modelsdf1)) {

    modelsdf1$clustermat[i, sample(x=dim(modelsdf1$clustermat)[2],
                                   size=initialclusters,
                                   replace=FALSE)] <- 1:initialclusters

  }
  colnames(modelsdf1$clustermat) <- unique(mal$placeid)
  # create sample adjacency
  adjacency2$cluster <- factor(fillbynearest(adjacency=adjacency, covariate=modelsdf1$clustermat[1,]))
  adjacency2$cluster <- NULL
  adjacency2$representative <- NULL

  # set up initial covariates
  envnames <- colnames(env)
  envnames <- envnames[!(envnames %in% c("placeid", "date", "year", "doy"))]
  envnames <- paste(envnames, "mat", sep="")

  x <- paste0("cov", seq(1:numofcovnts))

  for (j in x){
    # print(j)
    modelsdf1[j] <- sample(x=envnames,size=nrow(modelsdf1),
                           replace=TRUE)
  }

  modelsdf1$modelmeasure <- Inf
  modelsdf1$selectionprobability <- 1
  modelsdf1$generation <- 1

  # evaluate the first generation
  modelsdf1$modelmeasure <- evaluategeneration(models=modelsdf1,
                                               modeldata=mal,
                                               adjacency=adjacency,
                                               numofcovnts = numofcovnts)

  # rank models
  modelsdf1$selectionprobability <- rank(modelsdf1$modelmeasure, ties.method="random")
  modelsdf1$selectionprobability <- 1/modelsdf1$selectionprobability / sum(1/modelsdf1$selectionprobability, na.rm=TRUE)
  modelsdf1$selectionprobability[is.na(modelsdf1$selectionprobability)] <- 0

  # run the rest of the generations
  for (generation in 2:generations) {

    modelsdf1 <- addgeneration(models=modelsdf1,
                               modeldata=mal,
                               covariatenames=envnames,
                               individualspergeneration=individpergeneration,
                               adjacency=adjacency,
                               numofcovnts = numofcovnts)

    if (generation == generations) {

      write.csv(modelsdf1 %>% dplyr::select(!starts_with("clustermat", ignore.case = TRUE)),
                paste("generation_",
                      generation, ".csv",
                      sep=""))

    }


  }

}
