library(CINSignatureQuantification)

refitMixtureModel <- function(x = NULL,y = NULL,uninfPrior=TRUE){

    if( ncol(y) == 2 ) {
        # Poisson model
        x = round(x)
        if(uninfPrior) {
            postDatUnscaled = sapply(1:nrow(y), function(z) stats::dpois(x = x, lambda = y[[z,"Mean"]]) )
        } else {
            postDatUnscaled = sapply(1:nrow(y), function(z) stats::dpois(x = x, lambda = y[[z,"Mean"]]) *
                                         y[[z,"Weight"]])
        }

    } else {
        # Gaussian model
        if(uninfPrior) {
            postDatUnscaled = sapply(1:nrow(y), function(z) stats::dnorm(x = x, mean = y[[z,"Mean"]],
                                                                          sd = y[[z,"SD"]]) )
        } else {
            postDatUnscaled = sapply(1:nrow(y), function(z) stats::dnorm(x = x, mean = y[[z,"Mean"]],
                                                                          sd = y[[z,"SD"]]) * y[[z,"Weight"]] )
        }
    }
    return(postDatUnscaled)
}

convertSigs <- function(sigs = NULL,convertMat = NULL){
    if(is.null(sigs)){
        stop("requires signature matrix")
    }
    if(is.null(convertMat)){
        stop("requires signature conversion matrix")
    }
    # Match order of columns
    sigs = sigs[, match(rownames(convertMat), colnames(sigs)) ]
    # Convert sigs and rename components
    lifted = t( sigs %*% convertMat )
    return(lifted)
}

convertModelMap <- function(x = NULL,y = NULL,uninfPrior=TRUE,noCNFeat=TRUE){

    if(is.null(x)){
        stop("x is required - should be a list of mixture model components to be matched")
    }
    if(is.null(y)){
        stop("y is required - should be a list of mixture model components to match against")
    }

    stopifnot(is.logical(uninfPrior))
    stopifnot(is.logical(noCNFeat))

    allFeatures = names(y)

    if(noCNFeat){
        allFeatures <- allFeatures[!grepl(x = allFeatures,pattern = "copynumber")]
    }

    overlap <- lapply(allFeatures,FUN = function(feature){
        x_dat <- x[[feature]]$Mean
        y_dat <- y[[feature]]

        postDatUnscaled <- refitMixtureModel(x = x_dat,y = y_dat,uninfPrior = uninfPrior)
        postDatScaled <- data.frame(postDatUnscaled / rowSums(postDatUnscaled))

        colnames(postDatScaled) = paste0("Target_", feature, 1:ncol(postDatScaled))
        rownames(postDatScaled) <- paste0(feature, 1:nrow(postDatScaled))
        return(postDatScaled)
    })

    rownamesVec <- unlist(lapply(overlap,rownames))

    dtConv <- data.table::rbindlist(overlap, fill = TRUE)
    dtConv[ is.na(dtConv) ] <- 0
    mConv <- as.matrix(dtConv)
    rownames(mConv) <- rownamesVec

    return(mConv)
}

drewsMacintyreMapping <- function(){

    ## Load mixture models
    drewsMM = get(utils::data("Drews2022_TCGA_Mixture_Models",envir = environment()))
    macMM = get(utils::data("Macintyre2018_OV_Mixture_Models",envir = environment()))

    ## Load sig definitions
    drewsSigs <- get(utils::data("Drews2022_TCGA_Signatures",envir = environment()))
    macSigs <- get(utils::data("Macintyre2018_OV_Signatures",envir = environment()))

    ## Load Gold standard samples
    goldStandardData <- get(utils::data("TCGA_478_Samples_SNP6_GOLD",envir = environment()))

    ## Compute default acts
    drewsActs <- getActivities(quantifyCNSignatures(goldStandardData,method = "drews"))
    macActs <- getActivities(quantifyCNSignatures(goldStandardData,method = "mac"))

    ## remap drews signature defs to mac defs
    convertmap <- convertModelMap(x = drewsMM,y = macMM)
    convertMat <- convertSigs(sigs = drewsSigs,convertMat = convertmap)
    rownames(convertMat) <- gsub(x = rownames(convertMat),pattern = "Target_",replacement = "")

    ## Calculate cosine similarities
    actsim <- calculateCosineSim(x = drewsActs,macActs)
    sigsim <- calculateCosineSim(x = convertMat,y = macSigs[!grepl(pattern = "copynumber",rownames(macSigs)),])

    sigsimidx <- apply(sigsim,MARGIN = 2,FUN = function(z) ifelse(all(is.na(z)),NA,which(max(z,na.rm = T) == z)))
    mapTable <- data.frame()
    for(i in seq_len(length(sigsimidx))){
        idx <- sigsimidx[i]
        mapTable <- rbind(mapTable,c(names(sigsimidx)[i],rownames(sigsim)[idx],round(sigsim[idx,i],digits = 3)))
    }
    colnames(mapTable) <- c("drews","mac","cosine")

    sigsimidxActs <- apply(actsim,MARGIN = 2,FUN = function(z) ifelse(all(is.na(z)),NA,which(max(z,na.rm = T) == z)))
    mapTableActs <- data.frame()
    for(i in seq_len(length(sigsimidxActs))){
        idx <- sigsimidxActs[i]
        mapTableActs <- rbind(mapTableActs,c(names(sigsimidxActs)[i],rownames(actsim)[idx],round(actsim[idx,i],digits = 3)))
    }
    colnames(mapTableActs) <- c("drews","mac","cosine")

    sigsimidxinv <- apply(sigsim,MARGIN = 1,FUN = function(z) ifelse(all(is.na(z)),NA,which(max(z,na.rm = T) == z)))
    mapTableinv <- data.frame()
    for(i in seq_len(length(sigsimidxinv))){
        idx <- sigsimidxinv[i]
        mapTableinv <- rbind(mapTableinv,c(names(sigsimidxinv)[i],colnames(sigsim)[idx],round(sigsim[i,idx],digits = 3)))
    }
    colnames(mapTableinv) <- c("mac","drews","cosine")

    sigsimidxinvActs <- apply(sigsim,MARGIN = 1,FUN = function(z) ifelse(all(is.na(z)),NA,which(max(z,na.rm = T) == z)))
    mapTableinvActs <- data.frame()
    for(i in seq_len(length(sigsimidxinvActs))){
        idx <- sigsimidxinvActs[i]
        mapTableinvActs <- rbind(mapTableinvActs,c(names(sigsimidxinvActs)[i],colnames(actsim)[idx],round(actsim[i,idx],digits = 3)))
    }
    colnames(mapTableinvActs) <- c("mac","drews","cosine")

    return(list(map_drews_definitions=mapTable,
                map_mac_definitions=mapTableinv,
                cosine_definitions=sigsim,
                map_drews_activity=mapTableActs,
                map_mac_activity=mapTableinvActs,
                cosine_activity=actsim))
}

DrewsMacintyreMapping <- drewsMacintyreMapping()

rm(list=ls()[ls() != "DrewsMacintyreMapping"])
usethis::use_data(DrewsMacintyreMapping,overwrite = TRUE)
