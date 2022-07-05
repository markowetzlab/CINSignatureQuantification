## Create input matrix
calculateSampleByComponentMatrixDrews = function(brECNF, UNINFPRIOR = TRUE) {

    # Load mix models
    #allModels = get(load("data/Drews2022_TCGA_Mixture_Models.rda"))
    allModels = get(utils::data("Drews2022_TCGA_Mixture_Models",envir = environment()))
    allFeatures = names(allModels)

    # Loop over features and calculate posterior probabilities
    lMats = lapply(allFeatures, function(thisFeature) {

        thisEcnf = brECNF[[ thisFeature ]]
        thisModel = allModels[[ thisFeature ]]

        dat = as.numeric( thisEcnf[,2] )
        # We want a posterior, hence likelihood (density) times prior (weight)
        if( ncol(thisModel) == 2 ) {
            # Poisson model
            if(UNINFPRIOR){
                postDatUnscaled = sapply(1:nrow(thisModel), function(x) stats::dpois(x = dat, lambda = thisModel[[x,"Mean"]]) )
            } else {
                postDatUnscaled = sapply(1:nrow(thisModel), function(x) stats::dpois(x = dat, lambda = thisModel[[x,"Mean"]]) * thisModel[[x, "Weight"]] )
            }

        } else {
            # Gaussian model
            if(UNINFPRIOR){
                postDatUnscaled = sapply(1:nrow(thisModel), function(x) stats::dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]) )
            } else {
                postDatUnscaled = sapply(1:nrow(thisModel), function(x) stats::dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]) * thisModel[[x, "Weight"]] )
            }
        }

        # Normalise densities to probabilities
        postDatScaled = data.frame( postDatUnscaled / rowSums(postDatUnscaled) )
        postDatScaled$Sample = thisEcnf[,1]
        matSxC = stats::aggregate(. ~ Sample, postDatScaled, sum)
        rownames(matSxC) = matSxC$Sample
        matSxC$Sample = NULL
        matSxC = as.matrix(matSxC)

        # Should be sorted but just to be sure
        ## PS - single sample input results in a vector rather than a matrix
        ## This section reimplements as a matrix and renames the rows to match
        ## multi-sample input
        if(nrow(matSxC) == 1){
            s <- rownames(matSxC)
            matSxC <- t(as.matrix(matSxC[, order(thisModel$Mean) ]))
            rownames(matSxC) <- s
        } else {
            matSxC = matSxC[, order(thisModel$Mean) ]
        }
        colnames(matSxC) = paste0( thisFeature, 1:ncol(matSxC) )

        return(matSxC)

    } )

    # Return data
    allMats = do.call(cbind, lMats)
    lMats = list(sampleByComponent = allMats, model=allModels)
    return(lMats)
}
