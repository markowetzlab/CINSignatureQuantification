#' exportCIN
#'
#' This function exports data from a `CNQuant` or `SigQuant` class object.
#'
#' @param object CNQuant or SigQuant object
#' @param outputDir Output location provided as a writable file directory. This
#'   defaults to the current working directory.
#' @param outputPrefix A prefix added to the beginning of exported files.
#' @param includeExpName Include the object experiment name in the file name
#'   (default: FALSE).
#' @param sep Default file separator used in writing files (Default: '\\t')
#' @param fullExport Provide a full export of all data contained within the
#'   provided object, including reference data, feature values, models, and
#'   supporting information (default: FALSE).
#'
#' @return NULL
#' @export

exportCIN <- function(object,outputDir=NULL,outputPrefix=NULL,
                      includeExpName=FALSE,sep="\t",fullExport=FALSE){

    if(!any(methods::is(object) %in% c("CNQuant","SigQuant"))) {
        stop("input should be a CNQuant or SigQuant class object")
    }

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

    ## Include/exclude exp name
    if(includeExpName){
        expname <- object@ExpData@experimentName
    } else {
        expname <- ""
    }

    ## Set method name
    if(is.null(object@featFitting$method)){
        methodname <- ""
    } else {
        methodname <- object@featFitting$method
    }

    ## Save RDS
    rds.name <- fileNamer(outputDir,outputPrefix,expname,
                          methodname,fileName="CINQuantObj",
                          suffix=".rds")
    saveRDS(object = object,file = rds.name)

    ## Export segs
    writeData(table = getSegments(object),
        file = fileNamer(outputDir,outputPrefix,expname,methodname,
            fileName="cn_segments",suffix),
        sep = sep)

    writeData(table = rownamesToCol(getSamplefeatures(object)),
        file = fileNamer(outputDir,outputPrefix,expname,methodname,
            fileName="sample_info",suffix),
        sep = sep)

    ## Check and export features
    if(length(object@featData) > 0){
        writeData(table = exportFeats(getFeatures(object)),
            file = fileNamer(outputDir,outputPrefix,expname,methodname,
            fileName="cn_features",suffix),
            sep = sep)
    }

    if(length(object@featFitting) > 0){
        writeData(table = rownamesToCol(getSampleByComponent(object)),
            file = fileNamer(outputDir,outputPrefix,expname,methodname,
                fileName="cn_sampleByComponent",suffix),
            sep = sep)
    }

    ## check SigQuant class
    if(methods::is(object,"SigQuant")){
        writeData(table = rownamesToCol(getActivities(object)),
                  file = fileNamer(outputDir,outputPrefix,expname,methodname,
                                   fileName="sig_activities",suffix),
                  sep = sep)
    }

    if(fullExport){
        cat("full export does nothing currently\n")
    }
    cat(paste0("CIN data written to ",outputDir))
}

# Support function to write tabular files to disk
writeData <- function(table=NULL,file=NULL,sep="\t"){
    if(is.null(table)){
        stop("no data")
    }
    if(is.null(file)){
        stop("no file path")
    }

    data.table::fwrite(x = table,file = file,
                       append = FALSE,quote = FALSE,sep = sep,
                       row.names = FALSE,col.names = TRUE,na = NA)
}

# Support function to name files
fileNamer <- function(outputDir,outputPrefix,expname,methodname,fileName,suffix){
    outfile <- paste(Sys.Date(),outputPrefix,expname,methodname,fileName,sep = "_")
    outfile <- gsub(pattern = "__+",replacement = "_",x = outfile)
    outfile <- gsub(pattern = "^_",replacement = "",x = outfile)

    outputDir <- gsub(pattern = "//+",replacement = "/",paste0(outputDir,"/"))

    filePath <- paste0(outputDir,outfile,suffix)
    return(filePath)
}

# Converts CN feat list to tabular format
exportFeats <- function(x){
    ftab <- do.call(rbind,lapply(names(x)[1],FUN = function(y){
        feat <- x[[y]]
        feat$feature <- rep(y,times=nrow(feat))
        colnames(feat) <- c("ID","value","feature")
        return(feat)
    }))
    return(ftab)
}

# Retains rownames from matrix outputs for export to tabular format
rownamesToCol <- function(x) {
    x <- as.data.frame(x)
    rn <- rownames(x)
    rownames(x) <- NULL
    y <- cbind(rn,x)
    colnames(y) <- c("sample",colnames(x))
    return(y)
}
