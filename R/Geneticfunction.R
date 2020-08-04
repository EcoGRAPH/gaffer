


geneticimplement<- function(
  individpergeneration = NULL,
  initialclusters      = NULL,
  generations          =  NULL, 
  modeldata = mal,
  envdata = env,
  shapefile = shp,
  whichregion = whichregion,
  numofcovnts = numofcovnts
)
{
  
  adjacency <- shp
  
  
  if (!is.null(whichregion)) {
    
    adjacency <- adjacency[adjacency$R_NAME == whichregion,]  
    
  }
  
  adjacency2 <- adjacency
  adjacency <- adjacency[adjacency$NewPCODE %in% modeldata$placeid,]
  savenames <- adjacency$NewPCODE
  adjacency <- nb2listw(poly2nb(adjacency, queen=TRUE, row.names=adjacency$NewPCODE),
                        style="B")
  print(adjacency)
  adjacency2
  
  
# set up the first generation of models
modelsdf <- data.frame(generation=1,
                       mutation="firstgen",
                       modelnumber=1:individpergeneration,
                       modelmeasure=-Inf,
                       selectionprobability=1/individpergeneration)
print(modelsdf)
# set up initial clusters
modelsdf$clustermat <- matrix(rep(NA,
                                  size=nrow(modelsdf)*length(unique(mal$placeid)),
                                  replace=TRUE),
                              nrow=nrow(modelsdf),
                              ncol=length(unique(mal$placeid)))

modelsdf$clustermat[2,]
dim(modelsdf$clustermat)[2]
for (i in 1:nrow(modelsdf)) {
  
  modelsdf$clustermat[i, sample(x=dim(modelsdf$clustermat)[2],
                                size=initialclusters,
                                replace=FALSE)] <- 1:initialclusters
  
}
colnames(modelsdf$clustermat) <- unique(mal$placeid)
modelsdf$clustermat
unique(mal$placeid)
modelsdf
# create sample adjacency
adjacency2
adjacency2$cluster <- factor(fillbynearest(adjacency=adjacency, covariate=modelsdf$clustermat[1,]))

ggplot(adjacency2) + geom_sf(aes(fill=cluster))
adjacency2$representative <- modelsdf$clustermat[1,]
adjacency2$representative
ggplot(adjacency2) + geom_sf(aes(fill=cluster)) +
  theme(legend.position="none")
ggplot(adjacency2) + geom_sf(aes(fill=factor(representative))) +
  theme(legend.position="none")
adjacency2$cluster <- NULL
adjacency2$representative <- NULL

# set up initial covariates
envnames <- paste(envnames, "mat", sep="")

print(envnames)
modelsdf



numofcovnts<- 10


x <- paste0("cov", seq(1:numofcovnts))


for (j in x){
  print(j)
  modelsdf1[j] <- sample(x=envnames,size=nrow(modelsdf1),
                         replace=TRUE)
}

sample(x=envnames,size=nrow(modelsdf1),
       replace=TRUE)[1]
modelsdf1$cov1[1]

modelsdf1





modelsdf1$modelmeasure <- Inf
modelsdf1$selectionprobability <- 1
modelsdf1$generation <- 1
modelsdf
modelsdf1


# evaluate the first generation


modelsdf1$modelmeasure <- evaluategenerationf(models=modelsdf1,
                                              modeldata=mal,
                                              adjacency=adjacency,
                                              numberofcovnts = numofcovnts)
# rank models
modelsdf1$modelmeasure
modelsdf1$selectionprobability <- rank(modelsdf1$modelmeasure, ties.method="random")
modelsdf1$selectionprobability <- 1/modelsdf1$selectionprobability / sum(1/modelsdf1$selectionprobability, na.rm=TRUE)
modelsdf1$selectionprobability[is.na(modelsdf1$selectionprobability)] <- 0

# run the rest of the generations  
for (generation in 2:generations) {
  
  #modelsdf <- addgeneration(models=modelsdf,
  #                         modeldata=mal,
  #                          covariatenames=envnames,
  #                         individualspergeneration=individpergeneration,
  #                         adjacency=adjacency)
  modelsdf1 <- addgenerationf(models=modelsdf1,
                              modeldata=mal,
                              covariatenames=envnames,
                              individualspergeneration=individpergeneration,
                              adjacency=adjacency,
                              numberofcovnts = numofcovnts)
  
  
  if (generation %% 2 == 0) {
    
    write.csv(modelsdf1,
              paste(".\\outputs\\generation_",
                    generation, ".csv",
                    sep=""))
    
  }
  
}


# map best model
bestmodel <- modelsdf1[modelsdf1$modelmeasure == min(modelsdf1$modelmeasure, na.rm=TRUE),]
bestmodel <- bestmodel[1,]
bestmod <- data.frame(t(bestmodel$clustermat))
bestmod$NewPCODE <- rownames(bestmod)
names(bestmod) <- c("cluster", "NewPCODE")
bestmod$cluster <- factor(bestmod$cluster)

adjacency2$cluster <- bestmod$cluster
ggplot(adjacency2) + geom_sf(aes(fill=cluster))

adjacency2$fillcluster <- factor(fillbynearest(adjacency=adjacency,
                                               covariate=bestmodel[1,]$clustermat))
ggplot(adjacency2) + geom_sf(aes(fill=fillcluster))


}