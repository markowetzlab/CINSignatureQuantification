#' @rdname exportCIN-methods
#' @aliases exportCIN
setMethod("exportCIN",signature = "CNQuant",function(object,outputDir=NULL,
                                                     outputPrefix=NULL,
                                                     includeExpName=FALSE,
                                                     sep="\t",
                                                     fullExport=FALSE){
    cat("export stuff\n")
    ## Check outputDir
    if(!is.null(outputDir)){
        if(!dir.exists(outputDir)){
            stop("provided output directory does not exist")
        }
    } else {
        outputDir <- getwd()
    }

    ## Check prefix
    if(!is.null(outputPrefix)){
        outputPrefix <- as.character(outputPrefix)
    } else {
        outputPrefix <- ""
    }

    ## Check fullExport
    if(!is.logical(fullExport)){
        stop("fullExport argument should be logical - TRUE or FALSE")
    }

    ## default export
    ## Check CNQuant class
    if(sep == "\t"){
        suffix <- ".tsv"
    } else
    if(sep == ","){
        suffix <- ".csv"
    } else {
        suffix <- ".txt"
    }

    if(includeExpName){
        expname <- object@ExpData@experimentName
    } else {
        expname <- ""
    }

    if(inherits(object,what = "CNQuant")){
        cat("export segments...\n")
        writeData(table = getSegments(object),
                  file = fileNamer(outputDir,outputPrefix,expname,
                                   fileName="cn_segments",suffix),sep = sep)
        if(length(object@featData) > 0){
            cat("export features...\n")
            writeData(table = collapseFeats(getFeatures(object)),
                      file = fileNamer(outputDir,outputPrefix,expname,
                                       fileName="cn_features",suffix),sep = sep)
        }
    }

    ## check SigQuant class
    if(inherits(object,what = "SigQuant")){

    }

    if(fullExport){
        cat("full export")
    }
})

writeData <- function(table=NULL,file=NULL,sep="\t"){
    if(is.null(table)){
        stop("no data")
    }
    if(is.null(file)){
        stop("no file path")
    }

    print(file)

    # data.table::fwrite(x = table,file = file,
    #                    append = FALSE,quote = FALSE,sep = sep,row.names = FALSE,
    #                    col.names = TRUE,na = NA)
}

fileNamer <- function(outputDir,outputPrefix,expname,fileName,suffix){
    if(outputPrefix != "" & expname != ""){
        outfile <- paste(outputPrefix,expname,fileName,sep = "_")
    } else if(outputPrefix != "" & expname == ""){
        outfile <- paste(outputPrefix,fileName,sep = "_")
    } else if(outputPrefix == "" & expname != ""){
        outfile <- paste(expname,fileName,sep = "_")
    } else if(outputPrefix == "" & expname == ""){
        outfile <- fileName
    }

    outputDir <- gsub(pattern = "//+",replacement = "/",paste0(outputDir,"/"))

    filePath <- paste0(outputDir,outfile,suffix)
    return(filePath)
}

collapseFeats <- function(x){
    ftab <- do.call(rbind,lapply(names(x)[1],FUN = function(y){
        feat <- x[[y]]
        feat$feature <- rep(y,times=nrow(feat))
        colnames(feat) <- c("ID","value","feature")
        return(feat)
    }))
    return(ftab)
}

