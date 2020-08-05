
evaluategenerationf <- function(models=NULL,
                               modeldata=NULL,
                               adjacency=NULL,
                               numofcovnts = NULL) {
  
  
  # make sure we start from scratch
  models$modelmeasure <- Inf
  models$selectionprobability <- 1
  
  # create some helpful covariates
  modeldata$numdate <- as.numeric(modeldata$date)
  modeldata$doy     <- as.numeric(format(modeldata$date, "%j"))
  
  # set up the regression
  # myformula <- formula("objective ~ placeid + s(numdate, by=placeid, id=1) +
  #                       s(doy, bs='cc', id=2) +
  #                      cov1 + cov2 + cov3")
  # myfallbackformula <- formula("objective ~ s(numdate, id=1) +
  #                       s(doy, bs='cc', id=2) +
  #                      cov1 + cov2 + cov3")
  
  covts <- paste0("cov",seq(1:numofcovnts))
  
  
  #print(covts)
  
  z = as.factor(paste(covts,collapse="+"))
  
  y = paste("objective ~ placeid","s(numdate, by=placeid, id=1)","s(doy, id=2)",z, sep = "+")
  
 # print(y)
  myformula1 <- formula(y)
  
  #isTRUE(myformula == myformula1)
  
  
  #myformula <- formula("objective ~ placeid +
 #                       s(numdate, by=placeid, id=1) +
  #                      s(doy, id=2) +
  #                      cov1 + cov2 + cov3")
  
  
  myformula1
  
  
  y1 = paste("objective ~ s(numdate, id=1)","s(doy, id=2)",z, sep = "+")
  
  #print(y1)
  #myformula1 <- formula(y)
  myfallbackformula1 <- formula(y1)
    
    
  #myfallbackformula <- formula("objective ~ s(numdate, id=1) +
  #                     s(doy, id=2) +
  #                     cov1 + cov2 + cov3")
  

  myfallbackformula1
  
  # run batch_bam on all the models
  #print(models$clustermat)
  for (curmodelnum in 1:nrow(models)) {
    
    # put current clusters into model data
    curclusters <- data.frame(cluster = models$clustermat[curmodelnum,],
                              placeid = colnames(models$clustermat))
    #print(curclusters)
    curclusters$cluster <- fillbynearest(adjacency=adjacency,
                                         covariate=models$clustermat[curmodelnum,])
    modeldata <- left_join(modeldata,
                           curclusters,
                           by="placeid")
    
    # make sure our factors are indeed factors
    modeldata$cluster <- factor(modeldata$cluster)
    modeldata$placeid <- factor(modeldata$placeid)
    
    #print(modeldata)
    # figure out which covariates we need
    
    
    
    
   
   # print(models$cov1[curmodelnum])
   # print(head(models))
    for (j in covts){
     #print(models[,j][curmodelnum])
      modeldata[j] <- modeldata[,models[,j][curmodelnum]]
    }
    
    #print(head(modeldata))
    #modeldata$cov1 <- modeldata[,models$cov1[curmodelnum]]
    #modeldata$cov2 <- modeldata[,models$cov2[curmodelnum]]
    #modeldata$cov3 <- modeldata[,models$cov3[curmodelnum]]
    
    #print(head(modeldata))
    
    tryCatch({
      
      # fit the models
      modelfit <- batch_bam(data = modeldata,
                            bamargs = list("formula" = myformula1,
                                           "family" = gaussian(),
                                           "discrete" = TRUE,
                                           "nthread" = parallel::detectCores(logical=FALSE)-1),
                            bamargs_fallback = list("formula" = myfallbackformula1),
                            over = "cluster")
      #print(modelfit)
      myAICs <- extractAIC.batch_bam(models=modelfit)
      models$modelmeasure[curmodelnum] <- sum(myAICs[,2]) + (log(nrow(modeldata)) - 2)*sum(myAICs[,1])
      
      # models$modelmeasure[curmodelnum] <- sum(extractAIC.batch_bam(models=modelfit)[,2])
      # models$df <- sum(extract)
      
      # try cleaning up
      rm(modelfit)
      rm(curclusters)
      modeldata$cluster <- NULL
      gc()
      
    }, error = function(e) {
      
      models$modelmeasure[curmodelnum] <- Inf
      modeldata$cluster <- NULL
      rm(list=ls())
      gc()
      return(NULL)
      
    })
    
  }
  
  # retain smallest possible frame
  returnframe <- unlist(models$modelmeasure)
  rm(list=setdiff(ls(), "returnframe"))
  gc()
  
  return(returnframe)
  
}
