forciblyrevalue <- function(codes=NULL) {

  codes <- unlist(codes, use.names=FALSE)
  uniquevalues <- unique(codes[!is.na(codes)])
  newcodes <- match(x=codes,
                    table=uniquevalues,
                    nomatch=NA)
  return(newcodes)

}
