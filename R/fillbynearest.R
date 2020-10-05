fillbynearest <- function(shapefile=NULL,
                          covariate=NULL) {


  # get coordinates of centroids
  coords <- st_centroid(shapefile)
  mydistances <- st_distance(x=coords)

  # figure out which of the covariates are NA
  whichcovsnotNA <- which(!is.na(covariate))

  nearestnonNA <- rep(0, length(covariate))
  # for each row, figure out which non-NA node is nearest
  for (i in 1:length(nearestnonNA)) {

    nearestnonNA[i] <- which(mydistances[i,whichcovsnotNA] == min(mydistances[i,whichcovsnotNA]))[1]

  }
  nonNAvalues <- covariate[!is.na(covariate)]
  return(nonNAvalues[nearestnonNA])

}
