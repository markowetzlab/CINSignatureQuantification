#' addsampleFeatures
#'
#' `addsampleFeatures` adds custom sample-level data to the samplefeatData field
#' of a `CNQuant` or `SigQuant` class object. This can be additional sample
#' information (e.g. purity, tumour type, etc.) that can be stored within the
#' object and used in downstream analysis.
#'
#' @param object CNQuant or SigQuant class object
#' @param sample.data data.frame containing sample-level variables
#' @param id.col column containing sample identifiers (default: "sample")
#' @return CNQuant or SigQuant object with updated samplefeatData
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   data(test.sample.features)
#'   cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- addsampleFeatures(cnobj,sample.data=test.sample.features)
#'   getSamplefeatures(cnobj)
#' @export
addsampleFeatures <- function(object=NULL,sample.data=NULL,id.col = "sample"){
    if(is.null(object)){
        stop("No object provided")
    }
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
    if(!any(newDataSamples %in% rownames(sampFeat))){
        stop("no overlapping samples in sample.data")
    }
    mergedsampfeats <- merge.data.frame(sampFeat,sample.data,by.x = "row.names",
                                        by.y = id.col,all.x = T)
    rownames(mergedsampfeats) <- mergedsampfeats$Row.names
    mergedsampfeats <- mergedsampfeats[,-1]
    if(!all(rownames(mergedsampfeats) == names(object@segments))){
        stop("something terrible has happened with the data.frame order")
    }
    methods::initialize(object,samplefeatData=mergedsampfeats,
                        ExpData = methods::initialize(object@ExpData,
                                    last.modified = as.character(Sys.time())))
}
