extractCopynumberFeaturesDrews = function(CN_data, cores = 1, allowedError = 0.1, rmNorm = FALSE, build="hg19") {

    # Get chromosome length and centromere locations
    if (build == "hg19") {
        chrlen <- get(utils::data("hg19.chrom.sizes",envir = environment()))
        gaps <- get(utils::data("gap_hg19",envir = environment()))
    } else if (build == "hg38") {
        chrlen <- get(utils::data("hg38.chrom.sizes",envir = environment()))
        gaps <- get(utils::data("gap_hg38",envir = environment()))
    }
    centromeres = gaps[gaps[,8]=="centromere",]

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
                list(segsize = getSegsizeDrews(CN_data, rmNorm = rmNorm) )
            } else if (i == 2) {
                list(bp10MB = getBPnumDrews(CN_data,chrlen) )
            } else if (i == 3) {
                list(osCN = getOscillationDrews(CN_data,chrlen) )
            } else if (i == 4) {
                list(bpchrarm = getCentromereDistCountsDrews(CN_data,centromeres,chrlen) )
            } else if (i == 5) {
                list(changepoint = getChangepointCNDrews(CN_data, allowedError, rmNorm = rmNorm) )
            } else {
                # Technically not needed but kept for backwards compatibility of the code
                list(copynumber = getCNDrews(CN_data, rmNorm = rmNorm) )
            }
        }
        doParallel::stopImplicitCluster()

        # Another failsafe that the outcome is definitely numeric
        temp_list = unlist( temp_list, recursive = FALSE )
        outList = lapply(temp_list, function(thisDF) {
            thisDF[,2] = as.numeric(thisDF[,2])
            return(thisDF)
        })
        return( outList )

    } else {
        # Single core usage
        segsize<-getSegsizeDrews(CN_data,rmNorm = rmNorm)
        bp10MB<-getBPnumDrews(CN_data,chrlen)
        osCN<-getOscillationDrews(CN_data,chrlen)
        bpchrarm<-getCentromereDistCountsDrews(CN_data,centromeres,chrlen)
        changepoint<-getChangepointCNDrews(CN_data,rmNorm = rmNorm)
        copynumber<-getCNDrews(CN_data)

        temp_list = list(segsize=segsize,bp10MB=bp10MB,osCN=osCN,bpchrarm=bpchrarm,changepoint=changepoint,copynumber=copynumber)
        #temp_list = unlist( temp_list, recursive = FALSE )
        outList = lapply(temp_list, function(thisDF) {
            thisDF[,2] = as.numeric(thisDF[,2])
            return(thisDF)
        })
        return( outList )

    }

}
