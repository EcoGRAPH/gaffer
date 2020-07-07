rm(list=ls())

packages <- c("ggplot2", "plyr", "dplyr",
              "tidyverse", "reshape2", "pracma",
              "mgcv", "splines", "parallel", "splitstackshape",
              "data.table", "quantreg", "MASS",
              "sf", "maptools", "spdep", "igraph","gaffar")
for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, install.packages(package, repos = "http://cran.us.r-project.org"))
    library(package, character.only=T)
  }
}

library(clusterapply)

source("gaffer functions.R")
options(warn=-1)

# settings
# set up lag and regression information
laglen   <- 181
dlagdeg  <- 8
nprincomps <- 5

# load the harmonized weekly data
mal <- read.csv("robust weekly malaria.csv")
head(mal)

# make sure each observation corresponds to a time and a place
mal$date <- as.Date(mal$date, "%Y-%m-%d")
mal$doy  <- as.numeric(format(mal$date, "%j"))
mal$placeid <- mal$NewPCODE
mal$NewPCODE <- NULL
mal <- mal[!is.na(mal$placeid),]
head(mal)

# include the environmental data
env_brdf <- read.csv("Allbrdf2013-01-01to2019-12-31.csv", stringsAsFactors=TRUE)
env_prec <- read.csv("Allpresp2013-01-01to2019-12-31.csv", stringsAsFactors=TRUE)
env_surf <- read.csv("Allst2013-01-01to2019-12-31.csv", stringsAsFactors=TRUE)

# select a region
whichregion <- "Amhara"
woredanames <- env_brdf[c("NewPCODE",
                          "R_NAME",
                          "W_NAME",
                          "Z_NAME")]
table(env_brdf$R_NAME)
if (!is.null(whichregion)) {
  
  woredanames <- woredanames[woredanames$R_NAME == whichregion,]
  mal <- mal[mal$placeid %in% woredanames$NewPCODE,]
  
}

# create an adjacency matrix
adjacency <- sf::st_read("Eth_Admin_Woreda_2019_20200205.shp")
print(adjacency)
if (!is.null(whichregion)) {
  
  adjacency <- adjacency[adjacency$R_NAME == whichregion,]  
  
}
adjacency2 <- adjacency
adjacency <- adjacency[adjacency$NewPCODE %in% mal$placeid,]
savenames <- adjacency$NewPCODE
adjacency <- nb2listw(poly2nb(adjacency, queen=TRUE, row.names=adjacency$NewPCODE),
                      style="B")
print(adjacency)
# adjacencylist <- list()
# for (i in 1:length(adjacency$neighbours)) {
# 
#   templist <- adjacency$neighbours[[i]]
#   adjacencylist[[i]] <- unlist(templist)
#   
# }
# adjacencylist
         
env_brdf[c("X", "R_NAME", "W_NAME", "Z_NAME")] <- NULL
env_prec[c("X", "R_NAME", "W_NAME", "Z_NAME")] <- NULL
env_surf[c("X", "R_NAME", "W_NAME", "Z_NAME")] <- NULL
env <- left_join(env_brdf, env_prec, by=c("NewPCODE", "doy", "year"))
env <- left_join(env, env_surf, by=c("NewPCODE", "doy", "year"))
rm(env_brdf,
   env_prec,
   env_surf)
gc()
# create some variables
env$placeid <- env$NewPCODE
env$NewPCODE <- NULL
env$date <- as.Date(paste(env$year,
                          env$doy,
                          sep="-"),
                          "%Y-%j")

# make sure we have
env <- env[!is.na(env$placeid),]
env <- env[!is.na(env$date),]

# anomalize all environmental values
envnames <- colnames(env)
envnames <- envnames[!(envnames %in% c("placeid", "date", "year", "doy"))]
# for (curcol in envnames) {
# 
#   # get anomalies and fill missing
#   env$tempvar <- env[,curcol]
#   tempsum <- dplyr::summarize(group_by(env, placeid, doy),
#                               meantempvar = mean(tempvar, na.rm=TRUE))
#   env$tempvar <- NULL
#   env <- left_join(env, tempsum, by=c("placeid", "doy"))
#   env[,curcol] <- env[,curcol] - env$meantempvar
#   env[is.na(env[,curcol]),curcol] <- 0
#   env$meantempvar <- NULL
#   
#   # # rank  
#   # env[,curcol] <- rank(env[,curcol],
#   #                      ties.method="average",
#   #                      na.last="keep")
#   
# }
# env$year <- NULL
# env$doy  <- NULL

# # perform a princomp
# envprincomp <- princomp(env[,envnames],
#                         cor=TRUE,
#                         scores=TRUE,
#                         fix_sign=TRUE)
# plot(envprincomp)
# 
# # put scores back into environmental set
# env$princomp <- envprincomp$scores
# env[,envnames] <- NULL

# lag the environmental data
datalagger <- expand.grid(placeid=unique(mal$placeid),
                          date=unique(mal$date),
                          lag=seq(from=0, to=laglen-1, by=1))
datalagger$laggeddate <- datalagger$date-datalagger$lag
datalagger <- left_join(datalagger, env,
                        by=c("placeid"="placeid",
                             "laggeddate"="date"))

envnames
for (curenv in envnames) {

  # pivot environmental data by lag
  envlagged <- dcast(datalagger, placeid + date ~ lag, value.var=curenv)
  names(envlagged) <- paste(curenv,names(envlagged),sep="_")

  # join to malaria data
  mal <- left_join(x=mal, y=envlagged,
                   by=c("placeid"=paste(curenv,"placeid",sep="_"),
                        "date"=paste(curenv,"date",sep="_")))

  # turn these lagged data into matrices
  mal$tempvar <- as.matrix(mal[,grep(x=colnames(mal),
                                     pattern=paste(curenv, "_", sep=""),
                                     fixed=TRUE)])
  colnames(mal)[length(colnames(mal))] <- paste(curenv, "mat", sep="")
  
  # replace with summaries
  mal[,paste(curenv, "mat", sep="")] <- mal[,paste(curenv, "mat", sep="")] %*%
    as.matrix(bs(0:(laglen-1), intercept=TRUE))

  # and get rid of the lagged data columns
  mal[,grep(x=colnames(mal), pattern=paste(curenv, "_", sep=""), fixed=TRUE)] <- NULL

}
names(mal)

# # pivot and put back into main frame as a matrix
# for (curprincomp in 1:nprincomps) {
# 
#   # pivot environmental data by lag
#   datalagger$temp <- datalagger$princomp[,curprincomp]
#   princompdata <- dcast(datalagger, placeid + date ~ lag, value.var="temp")
#   names(princompdata) <- paste("princomp_",names(princompdata),sep="")
# 
#   # join to malaria data
#   mal <- left_join(x=mal, y=princompdata,
#                    by=c("placeid"="princomp_placeid",
#                         "date"="princomp_date"))
# 
#   # turn these lagged data into matrices
#   mal$tempvar <- as.matrix(mal[,grep(x=colnames(mal),
#                                                 pattern="princomp_",
#                                                 fixed=TRUE)])
#   colnames(mal)[length(colnames(mal))] <- paste("pc_", curprincomp, sep="")
#   
#   # and get rid of the lagged data columns
#   mal[,grep(x=colnames(mal), pattern="princomp_", fixed=TRUE)] <- NULL
# 
#   # summarize by basis
#   mal[,paste("pc_", curprincomp, sep="")] <- mal[,paste("pc_", curprincomp, sep="")] %*%
#     as.matrix(bs(0:(laglen-1), intercept=TRUE))
# 
# }
# names(mal)

# # create a lag matrix
# mal$lagmat <- mal$pc_1
# for (curlag in 0:(laglen-1)) {
# 
#   mal$lagmat[,curlag+1] <- curlag
# 
# }

# the genetic algorithm
# determine which variable we're modeling
mal$objective <- mal$robustified1

# set up genetic algorithm parameters
individpergeneration <- 25
initialclusters      <- 5
generations          <- 100

# set up the first generation of models
modelsdf <- data.frame(generation=1,
                       mutation="firstgen",
                       modelnumber=1:individpergeneration,
                       modelmeasure=-Inf,
                       selectionprobability=1/individpergeneration)

# set up initial clusters
modelsdf$clustermat <- matrix(rep(NA,
                                  size=nrow(modelsdf)*length(unique(mal$placeid)),
                                  replace=TRUE),
                              nrow=nrow(modelsdf),
                              ncol=length(unique(mal$placeid)))
for (i in 1:nrow(modelsdf)) {
  
  modelsdf$clustermat[i, sample(x=dim(modelsdf$clustermat)[2],
                                size=initialclusters,
                                replace=FALSE)] <- 1:initialclusters
  
}
colnames(modelsdf$clustermat) <- unique(mal$placeid)

# create sample adjacency
adjacency2$cluster <- factor(fillbynearest(adjacency=adjacency, covariate=modelsdf$clustermat[1,]))
ggplot(adjacency2) + geom_sf(aes(fill=cluster))
adjacency2$representative <- modelsdf$clustermat[1,]
ggplot(adjacency2) + geom_sf(aes(fill=cluster)) +
  theme(legend.position="none")
ggplot(adjacency2) + geom_sf(aes(fill=factor(representative))) +
  theme(legend.position="none")
adjacency2$cluster <- NULL
adjacency2$representative <- NULL

# set up initial covariates
envnames <- paste(envnames, "mat", sep="")
modelsdf$cov1 <- sample(x=envnames,
                        size=nrow(modelsdf),
                        replace=TRUE)
modelsdf$cov2 <- sample(x=envnames,
                        size=nrow(modelsdf),
                        replace=TRUE)
modelsdf$cov3 <- sample(x=envnames,
                        size=nrow(modelsdf),
                        replace=TRUE)

# set up the model measures
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
                            adjacency=adjacency)
  
  if (generation %% 2 == 0) {

    write.csv(modelsdf,
              paste(".\\outputs\\generation_",
                    generation, ".csv",
                    sep=""))

  }
  
}

# map best model
bestmodel <- modelsdf[modelsdf$modelmeasure == min(modelsdf$modelmeasure, na.rm=TRUE),]
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

# adjacency2$cluster
# adjacency2$fillcluster
# ggplot(adjacency2) + geom_point(aes(x=cluster, y=fillcluster))
