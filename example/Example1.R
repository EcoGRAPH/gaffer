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

options(warn=-1)

# settings
# set up lag and regression information
#laglen   <- 181
#dlagdeg  <- 8
#nprincomps <- 5

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
  print(woredanames)
  mal <- mal[mal$placeid %in% woredanames$NewPCODE,]
  
}
head(mal)
# create an adjacency matrix
shp <- sf::st_read("Eth_Admin_Woreda_2019_20200205.shp")


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
print(env)
# make sure we have
env <- env[!is.na(env$placeid),]
env <- env[!is.na(env$date),]



# Data lagging process
data <- dataprocessing(laglen   = 181,
                       dlagdeg  = 8,
                       nprincomps =5,
                       modeldata = mal,
                       env = env)


mal <- as.data.frame(data[1])
env <- as.data.frame(data[2])


mal$objective <- mal$robustified1
head(mal)

####Genetic algorithm calling

geneticimplement(
  individpergeneration = 10,
  initialclusters      = 4,
  generations          = 3,
  modeldata = mal,
  envdata = env,
  shapefile = shp,
  whichregion = "Amhara",
  numofcovnts = 3)
