rm(list=ls())

packages <- c("ggplot2", "plyr", "dplyr", "tidyr",
              "tidyverse", "reshape2", "pracma",
              "mgcv", "splines", "parallel", "splitstackshape",
              "data.table", "quantreg", "MASS",
              "sf", "maptools", "spdep", "igraph","gaffer",
              "clusterapply", "lwgeom", "dtwclust")
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

# decide which variable we're modeling
mal$objective <- mal$robustified2

# screen those which have very small counts
sumcases <- dplyr::summarise(group_by(mal, placeid),
                             sumobjective=sum(exp(objective), na.rm=TRUE))
lowcasethreshold <- quantile(sumcases$sumobjective, probs=c(0.15))
ggplot(sumcases) + geom_histogram(aes(x=sumobjective)) +
  geom_vline(xintercept=lowcasethreshold, linetype=2, color="red") +
  scale_x_log10()
sumcases <- sumcases[sumcases$sumobjective <= lowcasethreshold,]
mal <- mal[!(mal$placeid %in% sumcases$placeid),]

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
rm(woredanames)
# get the shapefile
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

# test for missing in environmental and mal matrices
sum(is.na(mal))
sum(is.na(env))

# save the names of the environmental data
envnames <- colnames(env)
rm(envframe)
gc()


# call the genetic algorithm
modelsdf <- geneticimplement(individpergeneration = 25,
                             initialclusters      = 10,
                             initialcovars        = 3,
                             generations          = 500,
                             modeldata = mal,
                             envnames = envnames,
                             shapefile = shp,
                             #forcecovariates="none",
                             slice = 5)
                             #restartfilename="C:\\home\\work\\davis\\gaffer\\csv outputs\\generation_345.csv")



# temporary debugging
modelsdf <- read.csv("C:\\home\\work\\davis\\gaffer\\csv outputs\\constrained environmentals 3\\generation_5.csv")

# load the outputs from the constrained set
constraineddf <- read.csv("C:\\home\\work\\davis\\gaffer\\csv outputs\\constrained environmentals 3\\generation_250.csv")

####### MODEL 1: BEST GAFFER MODEL #######
model1 <- modelsdf[which(modelsdf$modelmeasure == min(modelsdf$modelmeasure, na.rm=TRUE))[1],]
curclusterseeds <- unlist(strsplit(x=model1$clusterseeds,
                                   split=",",
                                   fixed=TRUE))
curclusters <- data.frame(placeid=rep("", length(curclusterseeds)/2),
                          cluster=rep(0 , length(curclusterseeds)/2))
for (i in 1:(length(curclusterseeds)/2)) {

  curclusters$placeid[i] <- curclusterseeds[2*i-1]
  curclusters$model1cluster[i] <- as.numeric(curclusterseeds[2*i])

}

# put this back into a map for filling
shp <- left_join(shp, curclusters, by="placeid")
shp$model1cluster <- factor(fillbynearest(shapefile=shp, covariate=shp$model1cluster))
ggplot(shp) + geom_sf(aes(fill=model1cluster))

# evaluate this model
mal <- left_join(mal, st_drop_geometry(shp[c("placeid", "model1cluster")]), by="placeid")
model1vars <- model1$covars[1]
model1vars <- unlist(strsplit(x=model1vars,
                              split=",",
                              fixed=TRUE))

# figure out which type of cyclicals it has
model1cyclicals <- model1$cyclicals[1]
if (model1cyclicals == "none") {

  model1baseformula <- "objective ~ placeid + s(numdate, by=placeid, bs='tp', id=1)"

}
if (model1cyclicals == "percluster") {

  model1baseformula <- "objective ~ placeid + s(numdate, by=placeid, bs='tp', id=1) + s(doy, bs='cc', id=2)"

}
if (model1cyclicals == "perplaceid") {

  model1baseformula <- "objective ~ placeid + s(numdate, by=placeid, bs='tp', id=1) + s(doy, bs='cc', by=placeid, id=2)"

}

# create the best model formula
model1formula <- paste("s(lagmat, by=",
                       model1vars,
                       "mat, bs='tp')",
                       sep="",
                       collapse="+")
model1formula <- as.formula(paste(model1baseformula, model1formula, sep="+"))

# create formulas in case those fail on some placeids
basefallback <- "objective ~ s(numdate, id=1) + s(doy, bs='cc', id=2)"
fallbackformula <- as.formula(basefallback)

# calculate the best fit from the gaffer model
mal$placeid <- factor(mal$placeid)
model1fit  <- batch_bam(data = mal,
                        bamargs = list("formula" = model1formula,
                                       "family" = gaussian(),
                                       "discrete" = TRUE,
                                       "nthread" = parallel::detectCores(logical=FALSE)-1),
                        bamargs_fallback = list("formula" = fallbackformula),
                        over = "model1cluster")






####### MODEL 2: BEST CONSTRAINED MODEL #######
model2 <- constraineddf[which(constraineddf$modelmeasure == min(constraineddf$modelmeasure, na.rm=TRUE))[1],]
curclusterseeds <- unlist(strsplit(x=model2$clusterseeds,
                                   split=",",
                                   fixed=TRUE))
curclusters <- data.frame(placeid=rep("", length(curclusterseeds)/2),
                          cluster=rep(0 , length(curclusterseeds)/2))
for (i in 1:(length(curclusterseeds)/2)) {

  curclusters$placeid[i] <- curclusterseeds[2*i-1]
  curclusters$model2cluster[i] <- as.numeric(curclusterseeds[2*i])

}

# put this back into a map for filling
shp <- left_join(shp, curclusters, by="placeid")
shp$model2cluster <- factor(fillbynearest(shapefile=shp, covariate=shp$model2cluster))
ggplot(shp) + geom_sf(aes(fill=model2cluster))

# evaluate this model
mal <- left_join(mal, st_drop_geometry(shp[c("placeid", "model2cluster")]), by="placeid")
model2vars <- model2$covars[1]
model2vars <- unlist(strsplit(x=model2vars,
                              split=",",
                              fixed=TRUE))

# figure out which type of cyclicals it has
model2cyclicals <- model2$cyclicals[1]
if (model2cyclicals == "none") {

  model2baseformula <- "objective ~ placeid + s(numdate, by=placeid, bs='tp', id=1)"

}
if (model2cyclicals == "percluster") {

  model2baseformula <- "objective ~ placeid + s(numdate, by=placeid, bs='tp', id=1) + s(doy, bs='cc', id=2)"

}
if (model2cyclicals == "perplaceid") {

  model2baseformula <- "objective ~ placeid + s(numdate, by=placeid, bs='tp', id=1) + s(doy, bs='cc', by=placeid, id=2)"

}

# create the best model formula
model2formula <- paste("s(lagmat, by=",
                       model2vars,
                       "mat, bs='tp')",
                       sep="",
                       collapse="+")
model2formula <- as.formula(paste(model2baseformula, model2formula, sep="+"))

# create formulas in case those fail on some placeids
basefallback <- "objective ~ s(numdate, id=1) + s(doy, bs='cc', id=2)"
fallbackformula <- as.formula(basefallback)

# calculate the best fit from the gaffer model
mal$placeid <- factor(mal$placeid)
model2fit  <- batch_bam(data = mal,
                        bamargs = list("formula" = model2formula,
                                       "family" = gaussian(),
                                       "discrete" = TRUE,
                                       "nthread" = parallel::detectCores(logical=FALSE)-1),
                        bamargs_fallback = list("formula" = fallbackformula),
                        over = "model2cluster")













####### MODEL 3: SINGLECLUSTER MODEL #######
mal$constantone <- factor(1)
model3fit <- batch_bam(data = mal,
                       bamargs = list("formula" = model2formula,
                                      "family" = gaussian(),
                                      "discrete" = TRUE,
                                      "nthread" = parallel::detectCores(logical=FALSE)-1),
                       bamargs_fallback = list("formula" = fallbackformula),
                       over = "constantone")









####### MODEL 4: SECULAR MODEL #######
model4formula  <- as.formula("objective ~ placeid + s(numdate, bs='tp', id=1) + s(doy, bs='cc', id=2)")
model4fit   <- batch_bam(data = mal,
                         bamargs = list("formula" = model4formula,
                                        "family" = gaussian(),
                                        "discrete" = TRUE,
                                        "nthread" = parallel::detectCores(logical=FALSE)-1),
                         bamargs_fallback = list("formula" = fallbackformula),
                         over = "model1cluster")










####### MODEL 5: K-MEANS MODEL #######
# remove the long-term trend from every placeid
# trendfit <- gam(objective ~ s(numdate, by=placeid), data=mal)
trendfit <- lm(objective ~ poly(numdate, 4)*placeid, data=mal)
mal$trendfit <- predict(trendfit, newdata=mal)
mal$trendres <- mal$objective - mal$trendfit

# calculate Fourier coefficients for the remaining stuff
fouriers <- data.frame()
for (curplaceid in unique(mal$placeid)) {

  thisplace <- mal[mal$placeid == curplaceid,]
  tempfft   <- data.frame(placeid=curplaceid,
                          fft=fft(thisplace[order(thisplace$numdate),]$trendres),
                          freq=1:nrow(thisplace))

  fouriers <- bind_rows(fouriers, tempfft)

}
fouriers$absfft <- abs(fouriers$fft)
ggplot(fouriers) + geom_line(aes(x=freq,
                                 y=absfft,
                                 group=placeid),
                             alpha=0.3)
naivefft <- fouriers[fouriers$freq %in% 1:10,]
naivefft <- pivot_wider(data=naivefft,
                        names_from=freq,
                        values_from=absfft,
                        names_prefix="fft_",
                        id_cols=placeid)
# keep copy of naivefft for later
naivefftcopy <- naivefft
# rescale the columns
for (colnum in 2:ncol(naivefft)) {

  naivefft[,colnum] <- scale(naivefft[,colnum])

}

wsses <- data.frame()
for (clusters in 1:20) {

  mykmeans <- kmeans(naivefft[,2:ncol(naivefft)],
                     clusters)

  wsses <- bind_rows(wsses,
                     data.frame(clusters=clusters,
                                wss=mykmeans$tot.withinss))

}

ggplot(wsses) + geom_point(aes(x=clusters, y=wss))
elbowed <-  kmeans(naivefft[,2:ncol(naivefft)], 13)
naivefft$model5cluster <- elbowed$cluster
naivefftcopy$model5cluster <- elbowed$cluster
shp <- left_join(shp, naivefft[c("placeid", "model5cluster")],
                 by="placeid")
shp$model5cluster <- factor(shp$model5cluster)

# reconstruct the time series
fouriers <- left_join(fouriers,
                      naivefftcopy[c("placeid", "model5cluster")],
                      by="placeid")
reconfft <- dplyr::summarise(group_by(fouriers, model5cluster, freq),
                             meanfft = mean(fft))

reconresiduals <- data.frame()
for (curcluster in unique(reconfft$model5cluster)) {

  thiscluster <- reconfft[reconfft$model5cluster == curcluster,]
  tempfft   <- data.frame(cluster=curcluster,
                          resid=fft(thiscluster$meanfft, inverse=TRUE))
  tempfft$nobs <- 1:nrow(tempfft)

  reconresiduals <- bind_rows(reconresiduals, tempfft)

}
# normalize
reconresiduals$resid <- Real(reconresiduals$resid) / max(reconresiduals$nobs)
reconresiduals$model5cluster <- factor(reconresiduals$cluster)

# add reconstructed residuals back to model
mal <- mal[order(mal$placeid, mal$numdate),]
mal$nobs <- rowid(mal$placeid)
mal <- left_join(mal, naivefftcopy[c("placeid", "model5cluster")],
                 by="placeid")
mal$model5cluster <- factor(mal$model5cluster)
mal <- left_join(mal, reconresiduals,
                 by=c("model5cluster"="model5cluster", "nobs"))

# reconstruct the observation
mal$model5pred <- mal$trendfit + mal$resid












####### MODEL 6: K-MEANS MODEL v2 #######

mycoefs <- data.frame()
for (curplaceid in unique(mal$placeid)) {

  tempdf <- mal[mal$placeid == curplaceid,]

  kmeansv2form <- formula("objective ~ s(numdate, bs='tp', id=1) + s(doy, bs='cc', id=2) + s(lagmat, by=lst_meanmat, bs='tp') + s(lagmat, by=ndwi6mat, bs='tp') + s(lagmat, by=totprecmat, bs='tp')")
  kmeansv2mod <- bam(formula=kmeansv2form,
                     discrete=TRUE,
                     data=tempdf,
                     nthreads=detectCores(logical=FALSE)-1)

  tempdf <- as.data.frame(coef(kmeansv2mod))
  tempdf$placeid <- curplaceid
  tempdf$variable <- rownames(tempdf)
  names(tempdf) <- c("value",
                     "placeid",
                     "variable")

  rownames(tempdf) <- NULL

  mycoefs <- bind_rows(mycoefs,
                       tempdf)

}
mycoefs <- mycoefs[mycoefs$variable != "(Intercept)",]
mycoefs <- mycoefs[!grepl(x=mycoefs$variable,
                          pattern="numdate",
                          fixed=TRUE),]

mycoefswide <- tidyr::pivot_wider(data=mycoefs,
                           names_from=variable,
                           values_from=value)
for (curcol in 2:ncol(mycoefswide)) {

  mycoefswide[,curcol] <- scale(mycoefswide[,curcol])

}
# perform a princomp
mycoefpc <- princomp(mycoefswide[,2:ncol(mycoefswide)],
                     cor=TRUE,
                     scores=TRUE)
mycoefswide_pc <- data.frame(mycoefpc$scores)
mycoefswide_pc$placeid <- mycoefswide$placeid

wsses <- data.frame()
for (clusters in 1:20) {

  mykmeans <- kmeans(mycoefswide_pc[,2:(ncol(mycoefswide_pc)-1)],
                     clusters)

  wsses <- bind_rows(wsses,
                     data.frame(clusters=clusters,
                                wss=mykmeans$tot.withinss))

}
ggplot(wsses) + geom_line(aes(x=clusters,
                              y=wss))

elbowed <- kmeans(mycoefswide[,2:ncol(mycoefswide_pc)], 6)
mycoefswide$model6cluster <- elbowed$cluster
shp <- left_join(shp, mycoefswide[c("placeid", "model6cluster")],
                 by="placeid")
mal <- left_join(mal, mycoefswide[c("placeid", "model6cluster")],
                 by="placeid")

# predict on full model
mal$placeid <- factor(mal$placeid)
mal$model6cluster <- factor(mal$model6cluster)
kmeansv2form <- formula("objective ~ placeid + s(numdate, by=placeid, bs='tp', id=1) + s(doy, bs='cc', id=2) + s(lagmat, by=lst_meanmat, bs='tp') + s(lagmat, by=ndwi6mat, bs='tp') + s(lagmat, by=totprecmat, bs='tp')")
kmeansv2formfall <- formula("objective ~ s(numdate, bs='tp', id=1) + s(doy, bs='cc', id=2) + s(lagmat, by=lst_meanmat, bs='tp') + s(lagmat, by=ndwi6mat, bs='tp') + s(lagmat, by=totprecmat, bs='tp')")
model6fit  <- batch_bam(data = mal,
                        bamargs = list("formula" = kmeansv2form,
                                       "family" = gaussian(),
                                       "discrete" = TRUE,
                                       "nthread" = parallel::detectCores(logical=FALSE)-1),
                        bamargs_fallback = list("formula" = kmeansv2formfall),
                        over = "model6cluster")





####### MODEL EVALUATION AND PLOTTING #######
distributedlags <- data.frame()
modelpreds      <- data.frame(placeid   = mal$place,
                              date      = mal$date,
                              objective = mal$objective)
modelgofs       <- data.frame()
shapefilegofs   <- NULL
# declare which models we're taking smooths from
whichmodels <- list("1_gaffer"=model1fit,
                    "2_constrained"=model2fit,
                    "3_singlecluster"=model3fit,
                    "4_secular"=model4fit,
                    "5_kmeans"=NA,
                    "6_kmeansv2"=model6fit)
# which overs are we using?
overs       <- list("1_gaffer"="model1cluster",
                    "2_constrained"="model2cluster",
                    "3_singlecluster"="constantone",
                    "4_secular"="model1cluster",
                    "5_kmeans"="model5cluster",
                    "6_kmeansv2"="model6cluster")

# figure out which models we shouldn't predict with batch_bam
# and which don't have distributed lags
donotpred <- list("5_kmeans"="model5pred")
noDLs     <- list("5_kmeans", "6_kmeansv2")





# figure out the predictions and compute goodness of fit statistics
for (whichmodel in names(whichmodels)) {

  if (!(whichmodel %in% names(donotpred))) {

    # load the model and its parameters
    thisfit <- whichmodels[[whichmodel]]
    thisover <- overs[[whichmodel]]

    # calculate the predictions for this model
    modelpreds[,paste(whichmodel, "pred", sep="_")] <- clusterapply::predict.batch_bam(models=thisfit,
                                                          predictargs=NULL,
                                                          over=thisover,
                                                          newdata=mal)

    # define a convenience
    modelpreds$temppred <- modelpreds[,paste(whichmodel, "pred", sep="_")]
    # get some universal goodness of fit statistics
    tempdf <- data.frame(model=whichmodel,
                         df =sum(extractAIC.batch_bam(thisfit)$X1),
                         AIC=sum(extractAIC.batch_bam(thisfit)$X2),
                         pearson=cor(modelpreds$objective,
                                     modelpreds$temppred,
                                     use="complete.obs",
                                     method="pearson"),
                         spearman=cor(modelpreds$objective,
                                      modelpreds$temppred,
                                      use="complete.obs",
                                      method="spearman"),
                         MAE_log=mean(abs(modelpreds$objective - modelpreds$temppred)),
                         MAE    =mean(abs(exp(modelpreds$objective) - exp(modelpreds$temppred))))

  } else {

    # obtain predictions from the model from somewhere else
    modelpreds[,paste(whichmodel, "pred", sep="_")] <- mal[,donotpred[[whichmodel]]]

    # define a convenience
    modelpreds$temppred <- modelpreds[,paste(whichmodel, "pred", sep="_")]

    # get some universal goodness of fit statistics
    tempdf <- data.frame(model=whichmodel,
                         df =NA,
                         AIC=NA,
                         pearson=cor(modelpreds$objective,
                                     modelpreds$temppred,
                                     use="complete.obs",
                                     method="pearson"),
                         spearman=cor(modelpreds$objective,
                                      modelpreds$temppred,
                                      use="complete.obs",
                                      method="spearman"),
                         MAE_log=mean(abs(modelpreds$objective - modelpreds$temppred)),
                         MAE    =mean(abs(exp(modelpreds$objective) - exp(modelpreds$temppred))))

  }

  modelgofs <- bind_rows(modelgofs, tempdf)

  # get fit statistics by placeid
  tempdf <- dplyr::summarise(group_by(modelpreds, placeid),
                             pearson=cor(objective, temppred,
                                         use="complete.obs",
                                         method="pearson"),
                             spearman=cor(objective, temppred,
                                          use="complete.obs",
                                          method="spearman"))
  tempdf$model <- whichmodel
  shapefilegofs <- bind_rows(shapefilegofs, tempdf)

  # delete convenience variable
  modelpreds$temppred <- NULL

  # figure out the distributed lags, if appropriate
  if (!(whichmodel %in% noDLs)) {

    # extract its distributed lags
    for (i in 1:length(thisfit)) {

      thisfit_plot <- plot.gam(thisfit[[i]], select=1)
      for (j in 1:length(thisfit_plot)) {

        if (grepl(x=thisfit_plot[[j]]$ylab,
                  pattern="lagmat",
                  fixed=TRUE)) {

            tempdf <- data.frame(model = whichmodel,
                                 cluster = names(thisfit)[i],
                                 variable = thisfit_plot[[j]]$ylab,
                                 x = thisfit_plot[[j]]$x,
                                 fit = thisfit_plot[[j]]$fit)

            distributedlags <- bind_rows(distributedlags, tempdf)

        }

      }

    }

  }

}
distributedlags$variable <- unlist(lapply(strsplit(x=distributedlags$variable,
                                            split=":",
                                            fixed=TRUE), "[[", 2))

# display distributed lags
distributedlags$model_cluster <- paste(distributedlags$model,
                                       distributedlags$cluster,
                                       sep="_")

ggplot(distributedlags) + geom_line(aes(x=x, y=fit,
                                        group=model_cluster,
                                        color=model),
                                    size=2) +
  facet_wrap(~variable, scales="free")

# display universal fits
modelpreds <- pivot_longer(data=modelpreds, cols=ends_with("_pred"))
ggplot(modelpreds) + geom_hex(aes(x=objective, y=value)) +
  geom_abline(slope=1, intercept=0, linetype=2, color="red") +
  facet_wrap(~name, scales="free")

# plot by model
reverseshp <- left_join(shapefilegofs,
                        shp,
                        by="placeid")
st_geometry(reverseshp) <- reverseshp$geometry
ggplot(reverseshp) + geom_sf(aes(fill=pearson)) +
  ggtitle("Pearson") +
  scale_fill_gradient2(low="darkblue",
                       mid="lightyellow",
                       high="darkred",
                       midpoint=0.5) +
  facet_wrap(~model)
ggplot(reverseshp) + geom_sf(aes(fill=spearman)) +
  ggtitle("Spearman") +
  scale_fill_gradient2(low="darkblue",
                       mid="lightyellow",
                       high="darkred",
                       midpoint=0.5) +
  facet_wrap(~model)

# # plot time series
# for (curplaceid in unique(modelpreds$placeid)) {
#
#   tempdf <- modelpreds[modelpreds$placeid == curplaceid,]
#   thisplot <- ggplot(tempdf) + geom_line(aes(x=date,
#                                              y=objective,
#                                              group=name),
#                                          color="black") +
#     geom_line(aes(x=date, y=value,
#                   group=name,
#                   color=name,
#                   linetype=name)) +
#     ggtitle(curplaceid) +
#     theme(legend.position="bottom")
#   ggsave(thisplot,
#          file=paste(".\\time series outputs\\",
#                curplaceid,
#                ".png", sep=""))
#
#
# }

# plot maps
stackshp <- data.frame()
tempovers <- overs[!is.na(overs)]
for (i in 1:length(tempovers)) {

  # get list of models
  tempdf <- mal[c("placeid", tempovers[[i]])]
  names(tempdf)[2] <- "cluster"

  # uniquify
  tempdf <- distinct(tempdf, .keep_all=TRUE)

  # add model to list
  tempdf$model <- names(tempovers)[i]
  stackshp <- bind_rows(stackshp, tempdf)

}
stackshp <- left_join(stackshp,
                      shp[c("placeid", "geometry")],
                      by="placeid")
st_geometry(stackshp) <- stackshp$geometry
stackshp$cluster <- factor(stackshp$cluster)
ggplot(stackshp) + geom_sf(aes(fill=cluster)) +
  facet_wrap(~model, ncol=2)

# plot time series
for (i in 1:length(overs)) {
<<<<<<< HEAD
<<<<<<< HEAD

  if (!is.na(overs[[i]])) {

      thisover <- overs[[i]]
      mal$toplot <- mal[,thisover]

      thisplot <- ggplot(mal) + geom_line(aes(x=date,
                                              y=objective-trendfit,
                                              group=placeid)) +
        facet_wrap(~toplot, ncol=2, scales="free") +
        ggtitle(names(overs)[i])
      plot(thisplot)

=======
=======
>>>>>>> e01a540cce625c4b9ec477cc480d608f22b247aa

  if (!is.na(overs[[i]])) {

      thisover <- overs[[i]]
      mal$toplot <- mal[,thisover]

      thisplot <- ggplot(mal) + geom_line(aes(x=date,
                                              y=objective-trendfit,
                                              group=placeid)) +
        facet_wrap(~toplot, ncol=2, scales="free") +
        ggtitle(names(overs)[i])
      plot(thisplot)

<<<<<<< HEAD
>>>>>>> parent of 8f94e0a... commit before RRSV
=======
>>>>>>> e01a540cce625c4b9ec477cc480d608f22b247aa
  }

}
