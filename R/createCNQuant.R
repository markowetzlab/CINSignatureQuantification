#' createCNQuant
#'
#' @param data Unrounded absolute copy number data
#' @param experimentName A name for the experiment (default: defaultExperiment)
#' @param build A genome build specified as either hg19 or hg38 (default: hg19)
#'
#' @return A CNQuant class object
#' @export createCNQuant
#'
createCNQuant <- function(data=NULL,experimentName = "defaultExperiment",build = "hg19"){
    if(is.null(data)){
        stop("no data provided\n")
    }
    supported_builds <- c("hg19","hg38")
    if(!build %in% supported_builds){
        stop(paste0("unknown build - supported builds: ",paste0(supported_builds,collapse = ", ")))
    }
    if(is.character(data)){
        if(!file.exists(data)){
            stop("File not found\n")
        }
        if(file.exists(data)){
            header <- colnames(data.table::fread(input = data,
                                                 header = T,
                                                 colClasses = c("character","numeric","numeric","numeric","character"),
                                                 nrows = 1))
            if(!any(header == c("chromosome","start","end","segVal","sample"))){
                stop("Header does not match the required naming")
            }
            segTable <- data.table::fread(input = data,
                                          header = T,
                                          colClasses = c("character","numeric","numeric","numeric","character"))
            if(checkSegValRounding(segTable$segVal)){
                warning("segVal appears to be rounded, copy number signatures require unrounded absolute copy numbers")
            }
            ## Binned inputs (fixed width not supported yet)
            # if(checkbinned(segTable)){
            #     #segTable <- getSegTable()
            #     #split(segTable,f = as.factor(segTable$sample))
            # } else {
            #     segTable <- split(segTable,f = as.factor(segTable$sample))
            # }
            ## Temp split until fixed bin input implemented
            segTable <- droplevels(segTable)
            segTable <- split(segTable,f = as.factor(segTable$sample))

            samplefeatData <- generateSampleFeatData(x = segTable)
            methods::new("CNQuant",segments = segTable,samplefeatData = samplefeatData,
                ExpData = methods::new("ExpQuant",
                              build = build,
                              samples.full = length(segTable),
                              samples.current = length(segTable),
                              experimentName = experimentName))
        }
    } else if("QDNAseqCopyNumbers" %in% class(data)){
        segTable <- getSegTable(x = data)
        if(checkSegValRounding(segTable$segVal)){
            warning("segVal appears to be rounded, copy number signatures require unrounded absolute copy numbers")
        }
        segTable <- droplevels(segTable)
        segTable <- split(segTable,f = as.factor(segTable$sample))
        samplefeatData <- generateSampleFeatData(x = segTable)
        methods::new("CNQuant",segments = segTable,samplefeatData = samplefeatData,
                     ExpData = methods::new("ExpQuant",
                          build = build,
                          samples.full = length(segTable),
                          samples.current = length(segTable),
                          experimentName = experimentName))
    } else if(is.data.frame(data)){
        header <- colnames(data)
        if(!all(header %in% c("chromosome","start","end","segVal","sample"))){
            stop("Header does not match the required naming ['chromosome','start','end','segVal','sample']")
        }
        segTable <- data
        if(checkSegValRounding(segTable$segVal)){
            warning("segVal appears to be rounded, copy number signatures were defined on unrounded absolute copy numbers, use caution when interpretting and comparing between rounded and unrounded inputs.")
        }
        ## Binned inputs (fixed width not supported yet)
        # if(checkbinned(segTable)){
        #     #segTable <- getSegTable()
        #     #split(segTable,f = as.factor(segTable$sample))
        # } else {
        #     segTable <- split(segTable,f = as.factor(segTable$sample))
        # }
        ## Temp split until fixed bin input implemented
        segTable <- droplevels(segTable)
        segTable <- split(segTable,f = as.factor(segTable$sample))
        samplefeatData <- generateSampleFeatData(x = segTable)
        methods::new("CNQuant",segments=segTable,samplefeatData = samplefeatData,
                     ExpData = methods::new("ExpQuant",
                          build = build,
                          samples.full = length(segTable),
                          samples.current = length(segTable),
                          experimentName = experimentName))
    } else {
        stop("Unknown input format\n")
    }
}
