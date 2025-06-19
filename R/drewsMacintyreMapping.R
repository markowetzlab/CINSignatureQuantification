refitMixtureModel <- function(x = NULL,y = NULL,uninfPrior=TRUE){

    if( ncol(y) == 2 ) {
        # Poisson model
        x = round(x)
        if(uninfPrior) {
            # Uninformative prior. As all weights would be the same, we can just drop it as we scale later anyways.
            postDatUnscaled = sapply(1:nrow(y), function(z) dpois(x = x, lambda = y[[z,"Mean"]]) )
        } else {
            postDatUnscaled = sapply(1:nrow(y), function(z) dpois(x = x, lambda = y[[z,"Mean"]]) *
                                         y[[z,"Weight"]])
        }

    } else {
        # Gaussian model
        if(uninfPrior) {
            # Uninformative prior. As all weights would be the same, we can just drop it as we scale later anyways.
            postDatUnscaled = sapply(1:nrow(y), function(z) dnorm(x = x, mean = y[[z,"Mean"]],
                                                                          sd = y[[z,"SD"]]) )
        } else {
            postDatUnscaled = sapply(1:nrow(y), function(z) dnorm(x = x, mean = y[[z,"Mean"]],
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
        stop("x is required - should be a list of mixture model components to match against")
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

## Test inverse mac --> drews
# drewsMM <- Drews2022_TCGA_Mixture_Models
# macMM <- Macintyre2018_OV_Mixture_Models
# drewsSigs <- Drews2022_TCGA_Signatures
# macSigs <- Macintyre2018_OV_Signatures
#
# convertMat <- convertSigs(macSigs,convertModelMap(x = macMM,y = drewsMM))
# rownames(convertMat) <- gsub(x = rownames(convertMat),pattern = "Target_",replacement = "")
#
# SIGobj.drews.mod <- SIGobj.drews
# SIGobj.drews.mod@backup.signatures <- t(convertMat)
#
# plotDefinitions(SIGobj.drews,normalise = T)
# plotDefinitions(SIGobj.drews.mod,normalise = T)
