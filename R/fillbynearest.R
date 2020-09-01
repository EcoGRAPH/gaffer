fillbynearest <- function(adjacency=NULL,
                          covariate=NULL) {

  # get list of shortest paths
  mygraph <- graph.adjacency(adjmatrix=listw2mat(adjacency),
                             mode="undirected")
  mypaths <- shortest.paths(graph=mygraph)

  # figure out which of the covariates are NA
  whichcovsNA <- which(is.na(covariate))

  # set path lengths to infinite to nodes which are NA
  mypaths[,whichcovsNA] <- Inf

  nearestnonNA <- rep(0, length(covariate))
  # for each row, figure out which non-NA node is nearest
  for (i in 1:length(nearestnonNA)) {

    nearestnonNA[i] <- which(mypaths[i,] == min(mypaths[i,]))[1]

  }
  return(unlist(covariate[nearestnonNA], use.names=FALSE))

}
