extractCopynumberFeaturesMac <- function(CN_data,cores = 1, build="hg19"){
    #chrlen <- get(load("data/hg19.chrom.sizes.rda"))
    if (build == "hg19") {
        chrlen <- get(data("hg19.chrom.sizes",envir = environment()))
        #gaps <- get(load("data/gap_hg19.rda"))
        gaps <- get(data("gap_hg19",envir = environment()))
    } else if (build == "hg38") {
        chrlen <- get(data("hg38.chrom.sizes",envir = environment()))
        #gaps <- get(load("data/gap_hg19.rda"))
        gaps <- get(data("gap_hg38",envir = environment()))
    }
    centromeres <- gaps[gaps[,8]=="centromere",]
    if(cores > 1) {
        if (!requireNamespace("doParallel", quietly = TRUE)) {
            stop(
                "Package \"doParallel\" must be installed to use multiple threads/cores.",
                call. = FALSE
            )
        }
        # Multi-core usage
        `%dopar%` <- foreach::`%dopar%`
        doParallel::registerDoParallel(cores)
        i <- NULL
        temp_list = foreach::foreach(i=1:6) %dopar% {
            if(i == 1){
                list(segsize = getSegsizeMac(CN_data) )
            } else if (i == 2) {
                list(bp10MB = getBPnumMac(CN_data,chrlen) )
            } else if (i == 3) {
                list(osCN = getOscillationMac(CN_data,chrlen) )
            } else if (i == 4) {
                list(bpchrarm = getCentromereDistCountsMac(CN_data,centromeres,chrlen) )
            } else if (i == 5) {
                list(changepoint = getChangepointCNMac(CN_data) )
            } else {
                list(copynumber = getCNMac(CN_data) )
            }

        }
        doParallel::stopImplicitCluster()
        unlist( temp_list, recursive = FALSE )
    } else {
        segsize <- getSegsizeMac(CN_data)
        bp10MB <- getBPnumMac(CN_data,chrlen)
        osCN <- getOscillationMac(CN_data,chrlen)
        bpchrarm <- getCentromereDistCountsMac(CN_data,centromeres,chrlen)
        changepoint <- getChangepointCNMac(CN_data)
        copynumber <- getCNMac(CN_data)

        list(segsize=segsize,
             bp10MB=bp10MB,
             osCN=osCN,
             bpchrarm=bpchrarm,
             changepoint=changepoint,
             copynumber=copynumber)
    }
}
