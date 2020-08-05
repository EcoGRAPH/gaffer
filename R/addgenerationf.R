forciblyrevalue <- function(codes=NULL) {
  
  newcodes <- codes
  
  if (is.data.frame(codes)) {
    
    tempcodes <- t(codes)
    for (curcode in 1:length(unique(sort(tempcodes)))) {
      
      toreplace <- unique(sort(tempcodes))[curcode]
      if (!is.na(toreplace)) {
        
        newcodes[codes == toreplace] <- curcode
        
      }
      
    }
    
  } else {
    
    for (curcode in 1:length(unique(sort(codes)))) {
      
      toreplace <- unique(sort(codes))[curcode]
      if (!is.na(toreplace)) {
        
        newcodes[codes == toreplace] <- curcode
        
      }
    }
    
  }
  return(newcodes)
  
}

addgenerationf <- function(models=NULL,
                          modeldata=NULL,
                          individualspergeneration=10,
                          covariatenames=NULL,
                          adjacency=NULL,
                          numofcovnts = NULL,
                          mutationprobabilities=c(0.60, 0.30, 0.05, 0.05),
                          mutationtypes=c("deletenode", "addnode", "recruit", "mutatevariables")) {
  
  # obtain the last generation
  lastgennum <- max(models$generation, na.rm=TRUE)
  lastgen <- models[models$generation == lastgennum,]
  
  # choose parents from the previous generation according to selectionprobability
  newgen <- lastgen[sample(x=nrow(lastgen),
                           size=individualspergeneration,
                           replace=TRUE,
                           prob=lastgen$selectionprobability),]
  
  # for each new child, select a possible mutation
  newgen$mutation <- sample(x=length(mutationprobabilities),
                            size=nrow(newgen),
                            replace=TRUE,
                            prob=mutationprobabilities)
  newgen$mutation <- mutationtypes[newgen$mutation]
  
  covts <- paste0("cov",seq(1:numofcovnts))
  
  
  # for each new child, mutate as chosen above
  for (i in 1:(nrow(newgen)-1)) {
    
    # if this is a merge, then pick two groups to merge
    if (newgen$mutation[i] == "deletenode") {
      
      maxgroup <- max(newgen$clustermat[i,], na.rm=TRUE)
      if (maxgroup >= 2) {
        
        whichnonNA <- which(!is.na(newgen$clustermat[i,]))
        whichtodelete <- whichnonNA[sample(x=length(whichnonNA),
                                           size=1)]
        newgen$clustermat[i, whichtodelete] <- NA
        
      } else {
        
        # if we can't perform a merge, then split instead
        newgen$mutation[i] <- "addnode"
        
      }
      
    }
    
    # if this is a split, then pick one group to split
    if (newgen$mutation[i] == "addnode") {
      
      maxgroup <- max(newgen$clustermat[i,], na.rm=TRUE)
      whichNA <- which(is.na(newgen$clustermat[i,]))
      whichtoadd <- whichNA[sample(x=length(whichNA), size=1)]
      newgen$clustermat[i, whichtoadd] <- (maxgroup+1)
      
    }
    
    # if this is a split, then pick one group to split
    if (newgen$mutation[i] == "recruit") {
      
      whichNA <- which(is.na(newgen$clustermat[i,]))
      whichnotNA <- which(!is.na(newgen$clustermat[i,]))
      whichtoadd <- whichNA[sample(x=length(whichNA), size=1)]
      whichbeingadded <- whichnotNA[sample(x=length(whichnotNA), size=1)]
      
      newgen$clustermat[i, whichtoadd] <- newgen$clustermat[i, whichbeingadded]
      
    }
    
    # if we're mutating variables
    if (newgen$mutation[i] == "mutatevariables") {
      
      # resample the covariates
      #print(newgen$cov1[i])
      for (j in covts){
       # print(j)
        z = newgen[,j]
        #print(z[i])
        z[i] <- sample(x=covariatenames, size=1)
      }
      
     # newgen$cov1[i] <- sample(x=covariatenames, size=1)
     # newgen$cov2[i] <- sample(x=covariatenames, size=1)
     # newgen$cov3[i] <- sample(x=covariatenames, size=1)
      
    } 
    
    # put back in order, make sure we're not missing anything
    newgen$clustermat[i,] <- forciblyrevalue(newgen$clustermat[i,])
    
  }
  
  # retain the alpha
  previousalpharow <- which(models$modelmeasure == min(models$modelmeasure, na.rm=TRUE))
  # make sure we're only retaining one
  previousalpharow <- previousalpharow[1]
  newgen[nrow(newgen),] <- models[previousalpharow,]
  
  # # get rid of the mutation signifier
  # newgen$mutation <- NULL
  
  # make sure we have model numbers
  newgen$modelnumber <- 1:nrow(newgen)
  
  # evaluate new generation
  newgen$modelmeasure <- evaluategenerationf(models=newgen,
                                            modeldata=mal,
                                            adjacency=adjacency,
                                            numofcovnts = numofcovnts)
  # rank models
  newgen$selectionprobability <- rank(newgen$modelmeasure, ties.method="random")
  newgen$selectionprobability <- 1/newgen$selectionprobability / sum(1/newgen$selectionprobability, na.rm=TRUE)
  newgen$selectionprobability[is.na(newgen$selectionprobability)] <- 0
  
  # update the generation counter
  newgen$generation <- lastgennum + 1
  
  # add this new generation to the population
  models <- rbind(models, newgen)
  
  # simplify rownames
  rownames(models) <- 1:nrow(models)
  
  return(models)
  
}
