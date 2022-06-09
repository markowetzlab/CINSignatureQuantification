setMethod("show", signature=c(object="CNQuant"),
          definition=function(object){
              cat(class(object)," object (initialised: ",object@ExpData@init.date,")\n\n",sep = "")
              cat("Experiment name:",object@ExpData@experimentName,"\n")
              cat("Segmented samples: ",object@ExpData@samples.current," (On init: ",object@ExpData@samples.full,")\n",sep = "")
              cat("CN feature data: ")
              if(length(object@featData) == 0){
                  cat("\n\tno data\n")
              } else {
                  cat("\n\tcn feature count: ",length(object@featData),"\n",sep = "")
                  cat("\tcn features: ",paste0(names(object@featData),collapse = ","),"\n",sep = "")
              }
              cat("Feature fitting: ")
              if(length(object@featFitting) == 0){
                  cat("\n\tno data\n")
              } else {
                  cat("\n\tsampleByComponent dim: ",dim(object@featFitting$sampleByComponent)[1]," x ",dim(object@featFitting$sampleByComponent)[2],"\n",sep = "")
                  cat("\tfitting method: ",object@featFitting$method,"\n",sep = "")
              }
              cat("Sample feature data:\n")
              cat("\tdim: ",dim(object@samplefeatData)[1]," x ",dim(object@samplefeatData)[2],"\n",sep = "")
              cat("\tfeatures: ",paste0(colnames(object@samplefeatData),collapse=","),"\n",sep = "")
              cat("Experiment:\n")
              cat("\tlast modified: ",object@ExpData@last.modified,"\n",sep = "")
              cat("\tgenome build: ",object@ExpData@build,"\n",sep = "")
              cat("`getExperiment(object)` for more details\n",sep = "")
})

#' Extract
#' @param x object to extract subset
#' @param i Indices of subset
#' @param j not used
#' @param ... not used
#' @param drop not used
#'
setMethod("[", signature=c("CNQuant", "numeric", "missing", "ANY"),
          definition=function(x, i, j, ..., drop=TRUE){
              segs <- x@segments[i]
              samplefeats <- x@samplefeatData[i,]
              if(length(x@featData) == 0){
                  feats <- x@featData
              } else {
                  subsetSamples <- names(segs)
                  feats <- subsetCNfeatures(x = x@featData,s = subsetSamples)
              }

              if(length(x@featFitting) == 0){
                  featsfit <- x@featFitting
              } else {
                  newfeatFitting <- subsetfeatFitting(x = x@featFitting,s = subsetSamples)
                  featsfit <- newfeatFitting
              }
              methods::initialize(x,
                         segments=segs,
                         featData=feats,
                         featFitting=featsfit,
                         samplefeatData=samplefeats,
                         ExpData = methods::initialize(x@ExpData,
                                              samples.current = length(segs),
                                              last.modified = as.character(Sys.time())))
          })

#' Extract
#' @param x object to extract subset
#' @param i Indices of subset
#' @param j not used
#' @param ... not used
#' @param drop not used
#'
setMethod("[", signature=c("CNQuant", "character", "missing", "ANY"),
          definition=function(x, i, j, ..., drop=TRUE){
              segs <- x@segments[names(x@segments) %in% i]
              samplefeats <- x@samplefeatData[rownames(x@samplefeatData) %in% i,]
              if(length(x@featData) == 0){
                  feats <- x@featData
              } else {
                  subsetSamples <- i
                  feats <- subsetCNfeatures(x = x@featData,s = subsetSamples)
              }
              if(length(x@featFitting) == 0){
                  featsfit <- x@featFitting
              } else {
                  newfeatFitting <- subsetfeatFitting(x = x@featFitting,s = subsetSamples)
                  featsfit <- newfeatFitting
              }
              methods::initialize(x,
                                  segments=segs,
                                  featData=feats,
                                  featFitting=featsfit,
                                  samplefeatData=samplefeats,
                                  ExpData = methods::initialize(x@ExpData,
                                                                samples.current = length(segs),
                                                                last.modified = as.character(Sys.time())))
          })
