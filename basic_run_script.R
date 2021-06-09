#
# Revised run script
# Modified by Dawn Nekorchuk
# last date : 2021-06-09

#reset
rm(list=ls())


## User settings --------------------------------------------------------------

#where to direct model final output (will be in datetime-stamped subfolder)
#note, while running, the individal gen output goes to the csv output and png output folders
output_base_folder <- file.path("model_outputs")


# select a region
whichregion <- "Amhara"


#settings for geneticimplement()
gaffer_settings <- c(individpergeneration = 20,  # how many individuals per genertaion? 20-ish is fine.
                     initialclusters      = 10,  # how many clusters do we start with? 10-ish is fine.
                     initialcovars        = 3,  # how many covariates do our first-generation models have?
                     generations          = 200,  # how many generations to run?
                     modeldata = mal,
                     envnames = envnames,
                     shapefile = shp,
                     forcecovariate="totprec,lst_day,ndwi6", # leave missing if we don't want to force covariates into each model
                     forcecyclicals="percluster",            # choose from percluster or perplaceid to force cyclicals
                     slice = 1)                              # how many generations between saves? 1 means save every generation. 5 means every fifth.
# if the GA crashed and needed to be restarted
# restartfilename = file.path("csv outputs", "generation_199.csv"))



## Package load ----------------------------------------------------------

#DMN: Note the following internal packages needed: "gaffer", "clusterapply"
devtools::install_github("EcoGRAPH/clusterapply@master")
#DMN uncomment if needed, but chances you may be developing and have a local build as this IS gaffer
#devtools::install_github("EcoGRAPH/gaffer@develop_ecograph")


packages <- c("rlang", "devtools", #DMN added
              "ggplot2", "plyr", "dplyr", "tidyr",
              "tidyverse", "reshape2", "pracma",
              "mgcv", "splines", "parallel", "splitstackshape",
              "data.table", "quantreg", "MASS",
              "sf", "maptools", "spdep", "igraph",
              "lwgeom", "dtwclust",
              "gaffer", "clusterapply") #Not CRAN, but should be installed from above, specifically listed last
for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, install.packages(package, repos = "http://cran.us.r-project.org"))
    library(package, character.only=T)
  }
}

options(warn=-1)


## Data import ----------------------------------------------------------------

# load the harmonized weekly data
mal_data_file <- file.path("data", "robust_weekly_malaria_w20200702.csv")
mal <- read.csv(mal_data_file)

# environmental data
env_brdf1 <- read_csv(file.path("data", "env_data/Export_Spectral_Data_2013-01-01_2013-12-31.csv"))
env_brdf2 <- read_csv(file.path("data", "env_data/Export_Spectral_Data_2014-01-01_2019-12-31.csv"))
env_brdf <- dplyr::bind_rows(env_brdf1, env_brdf2)

env_prec1 <- read_csv(file.path("data", "env_data/Export_Precip_Data_2013-01-01_2013-12-31.csv"))
env_prec2 <- read_csv(file.path("data", "env_data/Export_Precip_Data_2014-01-01_2019-12-31.csv"))
env_prec <- dplyr::bind_rows(env_prec1, env_prec2)

env_surf1 <- read_csv(file.path("data", "env_data/Export_LST_Data_2013-01-01_2013-12-31.csv"))
env_surf2 <- read_csv(file.path("data", "env_data/Export_LST_Data_2014-01-01_2019-12-31.csv"))
env_surf <- dplyr::bind_rows(env_surf1, env_surf2)

# shapefile
shp_file <- file.path("data", "Eth_Admin_Woreda_2019_20200702", "Eth_Admin_Woreda_2019_20200702.shp")
shp_full <- sf::st_read(shp_file)


## Data checks and setup ------------------------------------------------------

## Malaria Data

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


# pull woreda list that we have env data
woredanames <- env_brdf[c("NewPCODE",
                          "R_NAME",
                          "W_NAME",
                          "Z_NAME")]

# filter to region
if (!is.null(whichregion)) {

  woredanames <- woredanames[woredanames$R_NAME == whichregion,]
  mal <- mal[mal$placeid %in% woredanames$NewPCODE,]

}

# screen those which have very small counts
sumcases <- dplyr::summarise(group_by(mal, placeid),
                             sumobjective=sum(exp(objective), na.rm=TRUE))
# examine the cumulative distribution function of cases to find a pleasant cutoff
cdfframe <- data.frame(cases = seq(from=min(sumcases$sumobjective, na.rm=TRUE),
                                   to=1000,#max(sumcases$sumobjective, na.rm=TRUE),
                                   length.out=5000))
cdfframe$logcases <- log(cdfframe$cases + 1)
cdfframe$cdf <- 0
for (i in 1:nrow(cdfframe)) {

  cdfframe$cdf[i] <- mean(1*(sumcases$sumobjective <= cdfframe$cases[i]))

}
lowcasethreshold <- 600
ggplot(cdfframe) + geom_line(aes(x=cases,
                                 y=cdf)) +
  geom_vline(xintercept=lowcasethreshold, linetype=2, color="red") +
  scale_x_log10() +
  ggtitle("CDF of log total pfalc cases per woreda") +
  xlab("total pfalc cases")
ggplot(sumcases) + geom_histogram(aes(x=sumobjective), bins=20) +
  geom_vline(xintercept=lowcasethreshold, linetype=2, color="red") +
  scale_x_log10()

sumcases <- sumcases[sumcases$sumobjective <= lowcasethreshold,]

oldplaceids <- length(unique(mal$placeid))
mal <- mal[!(mal$placeid %in% sumcases$placeid),]
newplaceids <- length(unique(mal$placeid))

# what proportion of placeids have we retained?
newplaceids / oldplaceids
rm(woredanames)


## Shapefile

#set up
shp <- shp_full
shp$placeid <- shp$NewPCODE
shp$NewPCODE <- NULL
shp_region <- shp[shp$R_NAME == whichregion,]
#limits shp to only include woredas in the robustified data set
shp <- shp_region[shp_region$placeid %in% mal$placeid,]
plot(shp)


## Environmental

#join all env together
env_brdf[c("X", "R_NAME", "W_NAME", "Z_NAME")] <- NULL
env_prec[c("X", "R_NAME", "W_NAME", "Z_NAME")] <- NULL
env_surf[c("X", "R_NAME", "W_NAME", "Z_NAME")] <- NULL
env <- left_join(env_brdf, env_prec, by=c("NewPCODE", "doy", "year"))
env <- left_join(env, env_surf, by=c("NewPCODE", "doy", "year"))
rm(env_brdf,
   env_prec,
   env_surf)
gc()

# create some standard variables for use in the genetic algorithm
env$placeid <- env$NewPCODE
env$NewPCODE <- NULL
env$date <- as.Date(paste(env$year,
                          env$doy,
                          sep="-"),
                    "%Y-%j")

# make sure we have time and place
env <- env[!is.na(env$placeid),]
env <- env[!is.na(env$date),]

# only select environmental data from those districts which are going to be used
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


## Lag the environmental data ----------------------------------------------------------

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


## Genetic algorithm --------------------------------------------------------------------



#call the genetic algorithm
modelsdf <- geneticimplement(gaffer_settings)

# gaffer saves files in \\csv outputs\\ that contain information on cluster seeds, environmental
# covariates, model performance, etc. Any one of these csv files can be loaded instead of running
# the above function call.
# modelsdf <- read.csv("C:\\home\\work\\davis\\gaffer\\csv outputs\\21-02-01 - forced covariates on amhara per woreda\\generation_157.csv")


## Calculate additional results from geneticimplement return ---------------------------------------------


# recover some information about cluster size first, for debugging
modelsdf$maxcluster <- NA
for (currow in 1:nrow(modelsdf)) {

  theseseeds <- modelsdf$clusterseeds[currow]
  curclusterseeds <- unlist(strsplit(x=theseseeds,
                                     split=",",
                                     fixed=TRUE))
  curclusters <- data.frame(placeid=rep("", length(curclusterseeds)/2),
                            cluster=rep(0 , length(curclusterseeds)/2))
  for (i in 1:(length(curclusterseeds)/2)) {

    curclusters$placeid[i] <- curclusterseeds[2*i-1]
    curclusters$model1cluster[i] <- as.numeric(curclusterseeds[2*i])

  }
  modelsdf$maxcluster[currow] <- max(curclusters$model1cluster, na.rm=TRUE)

}

#DMN note: switching to more tidyverse style as I edit below

# adding ID field
models_tbl <- modelsdf %>%
  #tibble
  as_tibble() %>%
  #row # as ID
  mutate(modelID = row_number())


####### MODEL 1: BEST GAFFER MODEL #######
# obtain the model with the lowest AIC
# if multiple with same lowest value, takes the first (so, earliest generation)
#model1 <- modelsdf[which(modelsdf$modelmeasure == min(modelsdf$modelmeasure, na.rm=TRUE))[1],]
model1 <- models_tbl %>%
  slice_min(order_by = modelmeasure, n = 1, with_ties = FALSE)


# parse the string containing the cluster seeds so we know how many and where they are
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
shp_mod1 <- left_join(shp, curclusters, by="placeid")
shp_mod1$model1cluster <- factor(fillbynearest(shapefile=shp_mod1, covariate=shp_mod1$model1cluster))
shp_mod1$model1cluster <- factor(shp_mod1$model1cluster)
ggplot(shp_mod1) + geom_sf(aes(fill=model1cluster))


#plot with NA woredas too
p_mod1 <- ggplot() +
  geom_sf(data = shp_region) +
  geom_sf(data = shp_mod1, aes(fill = model1cluster)) +
  viridis::scale_fill_viridis(discrete = TRUE, option = "turbo") +
  theme_bw()

ggsave(plot = p_mod1, "cluster_map_model1.jpg",
       height = 6, width = 7, units = "in")

# evaluate this model
mal <- left_join(mal, st_drop_geometry(shp[c("placeid", "model1cluster")]), by="placeid")
model1vars <- model1$covars[1]
model1vars <- unlist(strsplit(x=model1vars,
                              split=",",
                              fixed=TRUE))

# figure out which type of cyclicals this model has
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



## Write out info and metadata/log file of best model ------------------------------------------------

# calc end date-time, to id results and prevent overwrites
res_dttm <- format(Sys.time(), "%Y%m%d_%H%M")

outpath <- file.path(output_base_folder, paste0(whichregion, "_", res_dttm))
dir.create(outpath, recursive = TRUE) #won't overwrite, if exists will do nothing except warn


# save model information for use in epidemia
code_outfile <- paste0("modelcodes_", res_dttm, ".csv")

write.csv(st_drop_geometry(shp_mod1),
          file = file.path(outpath, code_outfile))


#write log (here for now)
gs <- gaffer_settings[!grepl("^modeldata", names(gaffer_settings))]
gs <- gs[!grepl("^shapefile", names(gs))]

log <- list(gaffer_settings_trim = gs,
            whichregion = whichregion)

#save as rds for now
saveRDS(log,
        file = file.path(outpath, "log_parameters.RDS"))


#save models_tbl??
saveRDS(models_tbl,
        file = file.path(outpath, "models_tbl.RDS"))


#copies csv_outputs and model_outputs into final output folder
dir.create(file.path(outpath, "csv_outputs"), recursive = TRUE)
file.copy(from = list.files("csv outputs", full.names = TRUE),
          to = file.path(outpath, "csv_outputs"),
          overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = TRUE)

dir.create(file.path(outpath, "png_outputs"), recursive = TRUE)
file.copy(from = list.files("png outputs", full.names = TRUE),
          to = file.path(outpath, "png_outputs"),
          overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = TRUE)


## Evaluate and map several of the top performing models ----------------------------------

#get the top performing models
#this filters to anything within tolerance (tol = .Machine$double.eps^0.5) of the min score
mods_top <- filter(models_tbl, near(modelmeasure, min(models_tbl$modelmeasure)))

mods_top

mods_top$clusterseeds %>% unique()
#2021 06 09: all 31 have same clusterseeds, so same map





# # <<>> loop??
#
# # parse the string containing the cluster seeds so we know how many and where they are
# curclusterseeds <- unlist(strsplit(x=model1$clusterseeds,
#                                    split=",",
#                                    fixed=TRUE))
# curclusters <- data.frame(placeid=rep("", length(curclusterseeds)/2),
#                           cluster=rep(0 , length(curclusterseeds)/2))
# for (i in 1:(length(curclusterseeds)/2)) {
#
#   curclusters$placeid[i] <- curclusterseeds[2*i-1]
#   curclusters$model1cluster[i] <- as.numeric(curclusterseeds[2*i])
#
# }



# ####### MODEL EVALUATION AND PLOTTING #######
# distributedlags <- data.frame()
# modelpreds      <- data.frame(placeid   = mal$place,
#                               date      = mal$date,
#                               objective = mal$objective)
# modelgofs       <- data.frame()
# shapefilegofs   <- NULL
#
#
#
#
# # declare which models we're taking smooths from
# whichmodels <- list("model1"=model1fit,
#                     "model2"=model2fit)
# # which overs are we using?
# overs       <- list("model1"="model1cluster",
#                     "model2"="model2cluster")
#
# # figure out the predictions and compute goodness of fit statistics
# for (whichmodel in names(whichmodels)) {
#
#   # load the model and its parameters
#   thisfit <- whichmodels[[whichmodel]]
#   thisover <- overs[[whichmodel]]
#
#   # calculate the predictions for this model
#   modelpreds[,paste(whichmodel, "pred", sep="_")] <- clusterapply::predict.batch_bam(models=thisfit,
#                                                                                      predictargs=NULL,
#                                                                                      over=thisover,
#                                                                                      newdata=mal)
#
#   # define a convenience
#   modelpreds$temppred <- modelpreds[,paste(whichmodel, "pred", sep="_")]
#   # get some universal goodness of fit statistics
#   tempdf <- data.frame(model=whichmodel,
#                        df =sum(extractAIC.batch_bam(thisfit)$X1),
#                        AIC=sum(extractAIC.batch_bam(thisfit)$X2),
#                        pearson=cor(modelpreds$objective,
#                                    modelpreds$temppred,
#                                    use="complete.obs",
#                                    method="pearson"),
#                        spearman=cor(modelpreds$objective,
#                                     modelpreds$temppred,
#                                     use="complete.obs",
#                                     method="spearman"),
#                        totobs =sum(exp(modelpreds$objective)-1, na.rm=TRUE),
#                        MAE    =mean(abs(modelpreds$objective - modelpreds$temppred)))
#
#   modelgofs <- bind_rows(modelgofs, tempdf)
#
#   # get fit statistics by placeid
#   tempdf <- dplyr::summarise(group_by(modelpreds, placeid),
#                              pearson=cor(objective, temppred,
#                                          use="complete.obs",
#                                          method="pearson"),
#                              spearman=cor(objective, temppred,
#                                           use="complete.obs",
#                                           method="spearman"),
#                              MAE_log=mean(abs(objective-temppred), na.rm=TRUE),
#                              totobs=sum(exp(objective)-1, na.rm=TRUE))
#   tempdf$model <- whichmodel
#   shapefilegofs <- bind_rows(shapefilegofs, tempdf)
#
#   # delete convenience variable
#   modelpreds$temppred <- NULL
#
#   # extract its distributed lags
#   for (i in 1:length(thisfit)) {
#
#     thisfit_plot <- plot.gam(thisfit[[i]], select=0)
#     for (j in 1:length(thisfit_plot)) {
#
#       if (grepl(x=thisfit_plot[[j]]$ylab,
#                 pattern="lagmat",
#                 fixed=TRUE)) {
#
#         tempdf <- data.frame(model = whichmodel,
#                              cluster = names(thisfit)[i],
#                              variable = thisfit_plot[[j]]$ylab,
#                              x = thisfit_plot[[j]]$x,
#                              fit = thisfit_plot[[j]]$fit)
#
#         distributedlags <- bind_rows(distributedlags, tempdf)
#
#       }
#
#     }
#
#   }
#
# }
#
#
#
#
# distributedlags$variable <- unlist(lapply(strsplit(x=distributedlags$variable,
#                                                    split=":",
#                                                    fixed=TRUE), "[[", 2))
#
# # display distributed lags
# distributedlags$model_cluster <- paste(distributedlags$model,
#                                        distributedlags$cluster,
#                                        sep="_")
#
# ggplot(distributedlags) + geom_line(aes(x=x, y=fit,
#                                         group=model_cluster,
#                                         color=cluster,
#                                         linetype=model),
#                                     size=2) +
#   facet_wrap(~variable, scales="free")
#
#
#
#
# # display universal fits
# # ---------------------------Look into this-----------------------------------------------------
# modelpreds <- pivot_longer(data=modelpreds, cols=ends_with("_pred"))
# ggplot(modelpreds) + geom_hex(aes(x=objective, y=value)) +
#   geom_abline(slope=1, intercept=0, linetype=2, color="red") +
#   facet_wrap(~name, scales="free")
# # ----------------------------------------------------------------------------------------------
#
#
#
#
#
# # plot by model
# reverseshp <- left_join(shapefilegofs,
#                         shp,
#                         by="placeid")
# st_geometry(reverseshp) <- reverseshp$geometry
# ggplot() +
#   geom_sf(data=oldshp, fill="black") +
#   geom_sf(data=reverseshp, aes(fill=pearson)) +
#   ggtitle("Pearson") +
#   scale_fill_gradient2(low="darkblue",
#                        mid="lightyellow",
#                        high="darkred",
#                        midpoint=0.5) +
#   facet_wrap(~model) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.text=element_blank(),
#         axis.ticks=element_blank())
# ggplot(reverseshp) +
#   geom_sf(data=oldshp, fill="black") +
#   geom_sf(aes(fill=spearman)) +
#   ggtitle("Spearman") +
#   scale_fill_gradient2(low="darkblue",
#                        mid="lightyellow",
#                        high="darkred",
#                        midpoint=0.5) +
#   facet_wrap(~model) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.text=element_blank(),
#         axis.ticks=element_blank())
#
# # # plot time series
# # for (curplaceid in unique(modelpreds$placeid)) {
# #
# #   tempdf <- modelpreds[modelpreds$placeid == curplaceid,]
# #   thisplot <- ggplot(tempdf) + geom_line(aes(x=date,
# #                                              y=objective,
# #                                              group=name),
# #                                          color="black") +
# #     geom_line(aes(x=date, y=value,
# #                   group=name,
# #                   color=name,
# #                   linetype=name)) +
# #     ggtitle(curplaceid) +
# #     theme(legend.position="bottom")
# #   ggsave(thisplot,
# #          file=paste(".\\time series outputs\\",
# #                curplaceid,
# #                ".png", sep=""))
# #
# #
# # }
#
# # plot maps
# stackshp <- data.frame()
# tempovers <- overs[!is.na(overs)]
# for (i in 1:length(tempovers)) {
#
#   # get list of models
#   tempdf <- mal[c("placeid", tempovers[[i]])]
#   names(tempdf)[2] <- "cluster"
#
#   # uniquify
#   tempdf <- distinct(tempdf, .keep_all=TRUE)
#
#   # add model to list
#   tempdf$model <- names(tempovers)[i]
#   stackshp <- bind_rows(stackshp, tempdf)
#
# }
# stackshp <- left_join(stackshp,
#                       shp[c("placeid", "geometry")],
#                       by="placeid")
# st_geometry(stackshp) <- stackshp$geometry
# stackshp$cluster <- factor(stackshp$cluster)
# ggplot() + geom_sf(data=oldshp, fill="black") +
#   geom_sf(data=stackshp, aes(fill=cluster)) +
#   facet_wrap(~model, ncol=2) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.text=element_blank(),
#         axis.ticks=element_blank())
#
# # # plot time series
# # for (i in 1:length(overs)) {
# #
# #   if (!is.na(overs[[i]])) {
# #
# #       thisover <- overs[[i]]
# #       mal$toplot <- mal[,thisover]
# #
# #       thisplot <- ggplot(mal) + geom_line(aes(x=date,
# #                                               y=log(objective+1)-log(abs(trendfit)+1),
# #                                               group=placeid)) +
# #         facet_wrap(~toplot, ncol=2, scales="free") +
# #         ggtitle(names(overs)[i])
# #       plot(thisplot)
# #
# #   }
# #
# # }
#
# ggplot(shapefilegofs) + geom_point(aes(x=totobs, y=pearson)) +
#   scale_x_log10() +
#   facet_wrap(~model) +
#   ggtitle("Pearson rho per woreda vs. total pfalc malaria") +
#   xlab("total cases") + ylab("Pearson correlation")
# ggplot(shapefilegofs) + geom_point(aes(x=totobs, y=spearman)) +
#   scale_x_log10() +
#   facet_wrap(~model) +
#   ggtitle("Spearman rho per woreda vs. total pfalc malaria") +
#   xlab("total cases") + ylab("Spearman correlation")
# ggplot(shapefilegofs) + geom_point(aes(x=totobs, y=MAE_log)) +
#   scale_x_log10() +
#   facet_wrap(~model) +
#   ggtitle("MAE per woreda vs. total pfalc malaria") +
#   xlab("total cases") + ylab("MAE (log scale)")
