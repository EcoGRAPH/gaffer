addgeneration <- function(models=NULL,
                          modeldata=NULL,
                          individualspergeneration=10,
                          covariatenames=NULL,
                          adjacency=NULL,
                          placeids=NULL,
                          mutationprobabilities=c(0.40, 0.40, 0.10, 0.10),
                          mutationtypes=c("deletenode", "addnode", "dropvariable", "addvariable")) {

  # make absolutely sure all this is random
  set.seed(Sys.time())

  # obtain the last generation
  lastgennum <- max(models$generation, na.rm=TRUE)
  lastgen <- models[models$generation == lastgennum,]

  # choose parents from the previous generation according to selectionprobability
  newgen <- lastgen[sample(x=1:nrow(lastgen),
                           size=individualspergeneration,
                           replace=TRUE,
                           prob=lastgen$selectionprobability),]

  # for each new child, select a possible mutation
  newgen$mutation <- sample(x=1:length(mutationprobabilities),
                            size=nrow(newgen),
                            replace=TRUE,
                            prob=mutationprobabilities)
  newgen$mutation <- mutationtypes[newgen$mutation]

  # for each new child, mutate as chosen above
  for (i in 1:(nrow(newgen)-1)) {

    # figure out what the current model includes
    curcovars <- unlist(strsplit(x=newgen$covars[i],
                                 split=",",
                                 fixed=TRUE))

    curclusterseeds <- newgen$clusterseeds[i]
    splitseeds <- unlist(strsplit(x=curclusterseeds,
                         split=",",
                         fixed=TRUE))
    maxgroup <- length(splitseeds)/2

    # if this is a merge, then pick two groups to merge
    if (newgen$mutation[i] == "deletenode") {

      if (maxgroup >= 3) {

        # delete one of the seeds
        whichtodelete <- sample(x=1:(length(splitseeds)/2), size=1)
        splitseeds <- splitseeds[-c(2*whichtodelete-1, 2*whichtodelete)]

        # assign this new code
        newclusterseed <- paste(splitseeds, collapse=",")
        newgen$clusterseeds[i] <- newclusterseed

      } else {

        # if we can't perform a merge, then split instead
        newgen$mutation[i] <- "addnode"

      }

    }

    # if this is a split, then pick one group to split
    if (newgen$mutation[i] == "addnode") {

      # pick a new placeid that is not already present in our string
      whichadding <- placeids[which(!(placeids %in% splitseeds))]
      whichadding <- whichadding[sample(x=1:length(whichadding), size=1)]

      # add this to the list
      splitseeds <- c(splitseeds, whichadding, maxgroup+1)

      # assign this new code
      newclusterseed <- paste(splitseeds, collapse=",")
      newgen$clusterseeds[i] <- newclusterseed

    }

    # if we're deleting a variable
    if (newgen$mutation[i] == "dropvariable") {

      # figure out which variable we're deleting

      if (curcovars != "none") {

        whichtodelete <- sample(1:length(curcovars), size=1)
        curcovars <- curcovars[-c(whichtodelete)]

        if (length(curcovars) > 0) {

          # assign this new code
          newcovars <- paste(curcovars, collapse=",")
          newgen$covars[i] <- newcovars

        } else {

          newgen$covars[i] <- "none"

        }

      } else {

        newgen$covars[i] <- "none"

      }

    }

    # if we're adding a variable
    if (newgen$mutation[i] == "addvariable") {

      # figure out which variable we're adding
      whichtoadd <- covariatenames[!(covariatenames %in% curcovars)]

      if (length(whichtoadd) > 0) {

        whichtoadd <- whichtoadd[sample(1:length(whichtoadd), size=1)]

        if (curcovars != "none") {

          curcovars <- c(curcovars, whichtoadd)

        } else {

          curcovars <- whichtoadd

        }

        # assign these new covariates
        newcovars <- paste(curcovars, collapse=",")
        newgen$covars[i] <- newcovars

      } else {

        newgen$covars[i] <- curcovars

      }

    }

    # break apart once again to make sure the seeds are 1:n
    splitseeds <- unlist(strsplit(x=newgen$clusterseeds[i],
                                  split=",",
                                  fixed=TRUE))
    # get the cluster ids
    clusterids <- splitseeds[seq(from=2,
                                 to=length(splitseeds),
                                 by=2)]
    # revalue
    clusterids <- gaffer::forciblyrevalue(as.numeric(clusterids))

    # reassign
    splitseeds[seq(from=2,
                   to=length(splitseeds),
                   by=2)] <- as.character(clusterids)
    newgen$clusterseeds[i] <- paste(splitseeds,
                                    collapse=",")

  }

  # retain the alpha
  previousalpharow <- which(models$modelmeasure == min(models$modelmeasure, na.rm=TRUE))
  # make sure we're only retaining one
  previousalpharow <- previousalpharow[1]
  newgen[nrow(newgen),] <- models[previousalpharow,]
  newgen$mutation[nrow(newgen)] <- "alpha"

  # make sure we have model numbers
  newgen$modelnumber <- 1:nrow(newgen)

  # evaluate new generation
  newgen$modelmeasure <- evaluategeneration(models=newgen,
                                            modeldata=mal,
                                            adjacency=adjacency,
                                            placeids=placeids)
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
