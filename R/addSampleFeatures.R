#' addsampleFeatures
#'
#' Adds custom sample-level data to the samplefeatData field of a CNQuant or SigQuant object.
#' This can be additional sample information (purity, tumour type, etc.) that can
#' be used in downstream analysis.
#'
#' @param object CNQuant or SigQuant class object
#' @param sample.data data.frame containing sample-level variables
#' @param id.col column containing sample identifiers
#'
#' @return CNQuant or SigQuant object with updated samplefeatData
#' @export
#'
addsampleFeatures <- function(object,sample.data=NULL,id.col = "sample"){
    if(!class(object) %in% c("CNQuant","SigQuant")){
        stop("this function requires a CNQuant or SigQuant class object")
    }
    if(is.null(sample.data)){
        stop("no sample.data provided")
    }
    if(!is.data.frame(sample.data)){
        stop("sample data is not a data.frame")
    }
    if(!id.col %in% colnames(sample.data)){
        stop("id.col not found in sample.data")
    }
    sampFeat <- object@samplefeatData
    newDataSamples <- sample.data[,which(colnames(sample.data) == id.col)]
    if(!all(newDataSamples %in% rownames(sampFeat))){
        stop("no overlapping samples in sample.data")
    }
    mergedsampfeats <- merge.data.frame(sampFeat,sample.data,by.x = "row.names",by.y = id.col,all = T)
    rownames(mergedsampfeats) <- mergedsampfeats$Row.names
    mergedsampfeats <- mergedsampfeats[,-1]
    if(!all(rownames(mergedsampfeats) == names(object@segments))){
        stop("something terrible has happened with the data.frame order")
    }
    methods::initialize(object,samplefeatData=mergedsampfeats,
                        ExpData = methods::initialize(object@ExpData,
                                                      last.modified = as.character(Sys.time())))
}
