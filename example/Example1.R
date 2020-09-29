rm(list=ls())

packages <- c("ggplot2", "plyr", "dplyr",
              "tidyverse", "reshape2", "pracma",
              "mgcv", "splines", "parallel", "splitstackshape",
              "data.table", "quantreg", "MASS",
              "sf", "maptools", "spdep", "igraph","gaffer",
              "clusterapply")
for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, install.packages(package, repos = "http://cran.us.r-project.org"))
    library(package, character.only=T)
  }
}

options(warn=-1)

# load the harmonized weekly data
mal <- read.csv("robust weekly malaria.csv")

# make sure each observation corresponds to a time and a place
if (isDates(mal$date[1]) == FALSE) {
  mal$date <- as.Date(mal$date, format='%m/%d/%Y')
} else {
  mal$date <- as.Date(mal$date, format="%Y-%m-%d")
}

mal$doy  <- as.numeric(format(mal$date, "%j"))
mal$placeid <- mal$NewPCODE
mal$NewPCODE <- NULL
mal <- mal[!is.na(mal$placeid),]

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
if (!is.null(whichregion)) {

  woredanames <- woredanames[woredanames$R_NAME == whichregion,]
  mal <- mal[mal$placeid %in% woredanames$NewPCODE,]

}
# create an adjacency matrix
shp <- sf::st_read("Eth_Admin_Woreda_2019_20200205.shp")
shp$placeid <- shp$NewPCODE
shp$NewPCODE <- NULL
shp <- shp[shp$placeid %in% mal$placeid,]

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
# make sure we have time and place
env <- env[!is.na(env$placeid),]
env <- env[!is.na(env$date),]

# only select those which are going to be used
env <- env[env$placeid %in% mal$placeid,]

# create a frame that contains all the necessary environmental observations, missing if necessary
necessarydates <- expand.grid(date = unique(mal$date),
                              lag  = 0:181)
necessarydates$date <- necessarydates$date - necessarydates$lag
necessarydates$lag <- NULL
envframe <- expand.grid(placeid=unique(mal$placeid),
                        date   =unique(necessarydates$date))
env <- left_join(envframe, env, by=c("placeid", "date"))
rm(necessarydates)

# data lagging process
tempdf <- gaffer::dataprocessing(laglen = 181,
                                 modeldata = mal,
                                 env = env)
mal <- as.data.frame(tempdf[1])
env <- as.data.frame(tempdf[2])
rm(tempdf)

colnames(env)
colnames(mal)

# decide which variable we're modeling
mal$objective <- mal$robustified2

# call the genetic algorithm
modelsdf <- geneticimplement(individpergeneration = 2,
                             initialclusters      = 5,
                             initialcovars        = 1,
                             generations          = 4,
                             modeldata = mal,
                             envdata = env,
                             shapefile = shp,
                             slice = 2)#,
                             #restartfilename="C:\\home\\work\\davis\\gaffer\\csv outputs\\generation_10.csv")

# # load a saved file
# modelsdf <- read.csv("C:\\home\\work\\davis\\gaffer\\saved outputs\\20-09-18 - amhara pfalc anom\\generation_171.csv")

# reconstruct the best model
mybest <- modelsdf[which(modelsdf$modelmeasure == min(modelsdf$modelmeasure, na.rm=TRUE))[1],]
curclusterseeds <- unlist(strsplit(x=mybest$clusterseeds,
                                   split=",",
                                   fixed=TRUE))
curclusters <- data.frame(placeid=rep("", length(curclusterseeds)/2),
                          cluster=rep(0 , length(curclusterseeds)/2))
for (i in 1:(length(curclusterseeds)/2)) {

  curclusters$placeid[i] <- curclusterseeds[2*i-1]
  curclusters$cluster[i] <- as.numeric(curclusterseeds[2*i])

}

# put this back into a map
shp <- left_join(shp, curclusters, by="placeid")
adjacency <- shp[shp$placeid %in% unique(mal$placeid),]
placeids <- adjacency$placeids
adjacency <- nb2listw(poly2nb(adjacency,
                              queen=TRUE,
                              row.names=adjacency$placeid),
                      style="B")
shp$bestmodel <- factor(fillbynearest(adjacency=adjacency,
                                      covariate=shp$cluster))
ggplot(shp) + geom_sf(aes(fill=bestmodel))

# evaluate this model
bestmal <- left_join(mal, shp[c("placeid", "bestmodel")], by="placeid")
bestvars <- mybest$covars[1]
bestvars <- unlist(strsplit(x=bestvars,
                             split=",",
                             fixed=TRUE))
baseformula <- "objective ~ placeid + s(numdate, by=placeid, bs='tp', id=1) + s(doy, bs='cc', id=2)"
bestformula <- paste("s(",
                      bestvars,
                      "mat, by=lagmat, bs='tp')",
                      sep="",
                      collapse="+")
bestformula <- as.formula(paste(baseformula, bestformula, sep="+"))

basefallback <- "objective ~ s(numdate, id=1) + s(doy, bs='cc', id=2)"
fallbackformula <- as.formula(basefallback)

bestmal$placeid <- factor(bestmal$placeid)

modelfit <- batch_bam(data = bestmal,
                      bamargs = list("formula" = bestformula,
                                     "family" = gaussian(),
                                     "discrete" = TRUE,
                                     "nthread" = parallel::detectCores(logical=FALSE)-1),
                      bamargs_fallback = list("formula" = fallbackformula),
                      over = "bestmodel")

bestmal$bestpreds <- clusterapply::predict.batch_bam(models=modelfit,
                                                 predictargs=NULL,
                                                 over="bestmodel",
                                                 newdata=bestmal)
ggplot(bestmal) + geom_hex(aes(x=objective, y=bestpreds)) +
  geom_abline(slope=1, intercept=0, linetype=2, color="red") +
  ggtitle("preds vs. obs (best model) after 200 generations")

# calculate a null model fit
bestmal$constantone <- factor(1)
nullfit <- batch_bam(data = bestmal,
                     bamargs = list("formula" = bestformula,
                                    "family" = gaussian(),
                                    "discrete" = TRUE,
                                    "nthread" = parallel::detectCores(logical=FALSE)-1),
                     bamargs_fallback = list("formula" = fallbackformula),
                     over = "constantone")
# nullfit <- bam(data = bestmal,
#                formula = bestformula,
#                family = gaussian(),
#                discrete = TRUE,
#                nthread = parallel::detectCores(logical=FALSE)-1)

bestmal$nullpreds <- clusterapply::predict.batch_bam(models=nullfit,
                                                     predictargs=NULL,
                                                     over="constantone",
                                                     newdata=bestmal)

ggplot(bestmal) + geom_hex(aes(x=objective, y=nullpreds)) +
  geom_abline(slope=1, intercept=0, linetype=2, color="red") +
  ggtitle("preds vs. obs (null model) after 200 generations")

cor(x=bestmal$objective,
    y=bestmal$bestpreds,
    use="complete.obs")
cor(x=bestmal$objective,
    y=bestmal$nullpreds,
    use="complete.obs")

AIC1 <- sum(extractAIC.batch_bam(modelfit)$X2)
AIC0 <- sum(extractAIC.batch_bam(nullfit)$X2)

AIC1
AIC0
AIC1 - AIC0

# get summaries
bestsummary <- clusterapply::summary.batch_bam(modelfit)
nullsummary <- clusterapply::summary.batch_bam(nullfit)

# get estimated scales
bestscales <- unlist(lapply(bestsummary, "[[", "scale"))
nullscales <- unlist(lapply(nullsummary, "[[", "scale"))

bestscales
nullscales

exp(mean(log(bestscales)))
mean(bestscales)

# take a look at the summary
modelfit[[1]]

View(modelsdf)
