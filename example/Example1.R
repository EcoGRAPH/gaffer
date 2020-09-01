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

# data lagging process
tempdf <- gaffer::dataprocessing(laglen = 181,
                                 modeldata = mal,
                                 env = env)
mal <- as.data.frame(tempdf[1])
env <- as.data.frame(tempdf[2])
rm(tempdf)

# decide which variable we're modeling
mal$objective <- mal$robustified1

# call the genetic algorithm
modelsdf <- geneticimplement(individpergeneration = 30,
                             initialclusters      = 5,
                             initialcovars        = 1,
                             generations          = 200,
                             modeldata = mal,
                             envdata = env,
                             shapefile = shp)
