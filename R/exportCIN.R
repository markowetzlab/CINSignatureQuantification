#' @rdname exportCIN-methods
#' @aliases exportCIN
setMethod("exportCIN",signature = "CNQuant",function(object,outputDir=NULL,
                                                     outputPrefix=NULL,sep="\t",
                                                     fullExport=FALSE){
    cat("export stuff")
    if(!is.null(outputDir)){
        if(!dir.exists(outputDir)){
            stop("provided output directory does not exist")
        }
    } else {
        outputDir <- getwd()
    }

    if(!is.null(outputPrefix)){
        outputPrefix <- as.character(outputPrefix)
    } else {
        outputPrefix <- ""
    }
})
