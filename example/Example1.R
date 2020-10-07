rm(list=ls())

packages <- c("ggplot2", "plyr", "dplyr",
              "tidyverse", "reshape2", "pracma",
              "mgcv", "splines", "parallel", "splitstackshape",
              "data.table", "quantreg", "MASS",
              "sf", "maptools", "spdep", "igraph","gaffer",
              "clusterapply", "lwgeom")
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
lowcasethreshold <- quantile(sumcases$sumobjective, probs=c(0.10))
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

# # call the genetic algorithm
# modelsdf <- geneticimplement(individpergeneration = 2,
#                              initialclusters      = 5,
#                              initialcovars        = 1,
#                              generations          = 2,
#                              modeldata = mal,
#                              envdata = env,
#                              shapefile = shp,
#                              slice = 2)#,
#                              #restartfilename="C:\\home\\work\\davis\\gaffer\\csv outputs\\generation_10.csv")

# load a saved file
modelsdf <- read.csv("C:\\home\\work\\davis\\gaffer\\csv outputs\\generation_20.csv")

####### BEST GAFFER MODEL #######
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
shp$bestmodel <- factor(fillbynearest(shapefile=shp,
                                      covariate=shp$cluster))
ggplot(shp) + geom_sf(aes(fill=bestmodel))

# evaluate this model
bestmal <- left_join(mal, shp[c("placeid", "bestmodel")], by="placeid")
bestvars <- mybest$covars[1]
bestvars <- unlist(strsplit(x=bestvars,
                             split=",",
                             fixed=TRUE))

# figure out which type of cyclicals it has
bestcyclicals <- mybest$cyclicals[1]
if (bestcyclicals == "none") {

  baseformula <- "objective ~ placeid + s(numdate, by=placeid, bs='tp', id=1)"

}
if (bestcyclicals == "percluster") {

  baseformula <- "objective ~ placeid + s(numdate, by=placeid, bs='tp', id=1) + s(doy, bs='cc', id=2)"

}
if (bestcyclicals == "perplaceid") {

  baseformula <- "objective ~ placeid + s(numdate, by=placeid, bs='tp', id=1) + s(doy, bs='cc', by=placeid, id=2)"

}

# create the best model formula
bestformula <- paste("s(lagmat, by=",
                     bestvars,
                     "mat, bs='tp')",
                     sep="",
                     collapse="+")
bestformula <- as.formula(paste(baseformula, bestformula, sep="+"))

# create formulas in case those fail on some placeids
basefallback <- "objective ~ s(numdate, id=1) + s(doy, bs='cc', id=2)"
fallbackformula <- as.formula(basefallback)

# calculate the best fit from the gaffer model
bestmal$placeid <- factor(bestmal$placeid)
bestfit  <- batch_bam(data = bestmal,
                      bamargs = list("formula" = bestformula,
                                     "family" = gaussian(),
                                     "discrete" = TRUE,
                                     "nthread" = parallel::detectCores(logical=FALSE)-1),
                      bamargs_fallback = list("formula" = fallbackformula),
                      over = "bestmodel")

####### SINGLECLUSTER MODEL #######
bestmal$constantone <- factor(1)
singleclusterfit <- batch_bam(data = bestmal,
                              bamargs = list("formula" = bestformula,
                                             "family" = gaussian(),
                                             "discrete" = TRUE,
                                            "nthread" = parallel::detectCores(logical=FALSE)-1),
                              bamargs_fallback = list("formula" = fallbackformula),
                              over = "constantone")



# ####### SATURATED MODEL #######
# saturatedformula <- paste("s(",
#                           bestvars,
#                           "mat, by=lagmat, bs='tp')",
#                           sep="",
#                           collapse="+")
# saturatedformula  <- as.formula(paste("objective ~ s(numdate, bs='tp', id=1) + s(doy, bs='cc', id=2)",
#                                       saturatedformula, sep="+"))
# saturatedfit <- batch_bam(data = bestmal,
#                           bamargs = list("formula" = saturatedformula,
#                                          "family" = gaussian(),
#                                          "discrete" = TRUE,
#                                          "nthread" = parallel::detectCores(logical=FALSE)-1),
#                           bamargs_fallback = list("formula" = fallbackformula),
#                           over = "placeid")


####### SECULAR MODEL #######
secularformula  <- as.formula("objective ~ placeid + s(numdate, bs='tp', id=1) + s(doy, bs='cc', id=2)")
secularfit   <- batch_bam(data = bestmal,
                          bamargs = list("formula" = secularformula,
                                         "family" = gaussian(),
                                         "discrete" = TRUE,
                                         "nthread" = parallel::detectCores(logical=FALSE)-1),
                          bamargs_fallback = list("formula" = fallbackformula),
                          over = "bestmodel")

####### K-MEANS MODEL #######
# remove the long-term trend from every placeid
trendfit <- lm(objective ~ poly(numdate, 4)*placeid,
               data=bestmal)
bestmal$trendfit <- predict(trendfit, newdata=bestmal)
bestmal$trendres <- bestmal$objective - bestmal$trendfit
ggplot(bestmal) + geom_hex(aes(x=date, y=trendres)) +
  ggtitle("residual after removing quartic per woreda")
# calculate Fourier coefficients for the remaining stuff
fouriers <- data.frame()
for (curplaceid in unique(bestmal$placeid)) {

  thisplace <- bestmal[bestmal$placeid == curplaceid,]
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
elbowed <-  kmeans(naivefft[,2:ncol(naivefft)],
                   13)
naivefft$fftcluster <- elbowed$cluster
naivefftcopy$fftcluster <- elbowed$cluster
shp <- left_join(shp, naivefft[c("placeid", "fftcluster")],
                 by="placeid")
shp$fftcluster <- factor(shp$fftcluster)
ggplot(shp) + geom_sf(aes(fill=fftcluster))

# reconstruct the time series
fouriers <- left_join(fouriers,
                      naivefftcopy[c("placeid", "fftcluster")],
                      by="placeid")
reconfft <- dplyr::summarise(group_by(fouriers, fftcluster, freq),
                             meanfft = mean(fft))
reconresiduals <- data.frame()
for (curcluster in unique(reconfft$fftcluster)) {

  thiscluster <- reconfft[reconfft$fftcluster == curcluster,]
  tempfft   <- data.frame(cluster=curcluster,
                          resid=fft(thiscluster$meanfft, inverse=TRUE))
  tempfft$nobs <- 1:nrow(tempfft)

  reconresiduals <- bind_rows(reconresiduals, tempfft)

}
reconresiduals$resid <- Real(reconresiduals$resid)
reconresiduals$cluster <- factor(reconresiduals$cluster)

ggplot(reconresiduals) + geom_line(aes(x=nobs, y=resid, group=cluster, color=cluster))

# add reconstructed residuals back to model
bestmal <- bestmal[order(bestmal$placeid, bestmal$numdate),]
bestmal$nobs <- rowid(bestmal$placeid)
bestmal <- left_join(bestmal, naivefftcopy[c("placeid", "fftcluster")],
                     by="placeid")
bestmal$fftcluster <- factor(bestmal$fftcluster)
bestmal <- left_join(bestmal, reconresiduals,
                     by=c("fftcluster"="cluster", "nobs"))

# reconstruct the observation
bestmal$kmeans_pred <- bestmal$trendfit + bestmal$resid

ggplot(bestmal) + geom_hex(aes(x=objective, y=kmeans_pred)) +
  geom_abline(slope=1, intercept=0, linetype=2, color="red")

# # get summaries
# bestsummary <- clusterapply::summary.batch_bam(modelfit)
# singleclustersummary <- clusterapply::summary.batch_bam(singleclusterfit)
#
# # get estimated scales
# bestscales <- unlist(lapply(bestsummary, "[[", "scale"))
# singleclusterscales <- unlist(lapply(singleclustersummary, "[[", "scale"))

# # look at correlations by woreda
# gofbyplace <- dplyr::summarise(group_by(bestmal, placeid),
#                                bestcor = cor(objective, bestpreds, use="complete.obs"),
#                                singleclustercor = cor(objective, singleclusterpreds, use="complete.obs"))
# shp <- left_join(shp, gofbyplace, by="placeid")
# ggplot(shp) + geom_sf(aes(fill=bestcor)) +
#   scale_fill_gradient2(low="darkblue", mid="lightyellow", high="darkred", midpoint=0.5)
# ggplot(shp) + geom_sf(aes(fill=singleclustercor)) +
#   scale_fill_gradient2(low="darkblue", mid="lightyellow", high="darkred", midpoint=0.5)
#
# ggplot(shp) + geom_point(aes(x=singleclustercor, y=bestcor)) +
#   geom_abline(slope=1, intercept=0, linetype=2, color="red") +
#   xlim(0,1) + ylim(0,1) +
#   xlab("singlecluster performance") + ylab("gaffer performance")

####### MODEL EVALUATION AND PLOTTING #######
distributedlags <- data.frame()
modelpreds      <- data.frame(placeid   = bestmal$place,
                              date      = bestmal$date,
                              objective = bestmal$objective)
modelgofs       <- data.frame()
shapefilegofs   <- NULL
# declare which models we're taking smooths from
whichmodels <- list("singlecluster"=singleclusterfit,
                    "secular"=secularfit,
                    "gaffer"=bestfit,
                    "kmeans"=kmeans)
# this is really ugly - the batch_bam object should store its over
overs       <- list("gaffer"="bestmodel",
                    "singlecluster"="constantone",
                    "secular"="bestmodel",
                    "kmeans"=nullspace())
for (whichmodel in names(whichmodels)) {

  if (whichmodel != "kmeans") {

    # load the model and its parameters
    thisfit <- whichmodels[[whichmodel]]
    thisover <- overs[[whichmodel]]

    # calculate the predictions for this model
    modelpreds[,paste(whichmodel, "pred", sep="_")] <- clusterapply::predict.batch_bam(models=thisfit,
                                                          predictargs=NULL,
                                                          over=thisover,
                                                          newdata=bestmal)

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
    modelgofs <- bind_rows(modelgofs, tempdf)

  } else {

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
    modelgofs <- bind_rows(modelgofs, tempdf)

  }

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

  if (whichmodel != "kmeans") {

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

# plot time series
for (curplaceid in unique(modelpreds$placeid)) {

  head(modelpreds)
  tempdf <- modelpreds[modelpreds$placeid == curplaceid,]
  thisplot <- ggplot(tempdf) + geom_line(aes(x=date,
                                             y=objective,
                                             group=name),
                                         color="black") +
    geom_line(aes(x=date, y=value,
                  group=name,
                  color=name,
                  linetype=name)) +
    ggtitle(curplaceid) +
    theme(legend.position="bottom")
  ggsave(thisplot,
         file=paste(".\\time series outputs\\",
               curplaceid,
               ".png", sep=""))


}

View(modelgofs)
