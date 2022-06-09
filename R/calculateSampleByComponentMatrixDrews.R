## Create input matrix
calculateSampleByComponentMatrixDrews = function(brECNF, UNINFPRIOR = TRUE) {

    # Load mix models
    #allModels = get(load("data/Drews2022_TCGA_Mixture_Models.rda"))
    allModels = get(data("Drews2022_TCGA_Mixture_Models",envir = environment()))
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
                postDatUnscaled = sapply(1:nrow(thisModel), function(x) dpois(x = dat, lambda = thisModel[[x,"Mean"]]) )
            } else {
                postDatUnscaled = sapply(1:nrow(thisModel), function(x) dpois(x = dat, lambda = thisModel[[x,"Mean"]]) * thisModel[[x, "Weight"]] )
            }

        } else {
            # Gaussian model
            if(UNINFPRIOR){
                postDatUnscaled = sapply(1:nrow(thisModel), function(x) dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]) )
            } else {
                postDatUnscaled = sapply(1:nrow(thisModel), function(x) dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]) * thisModel[[x, "Weight"]] )
            }
        }

        # Normalise densities to probabilities
        postDatScaled = data.frame( postDatUnscaled / rowSums(postDatUnscaled) )
        postDatScaled$Sample = thisEcnf[,1]
        matSxC = aggregate(. ~ Sample, postDatScaled, sum)
        rownames(matSxC) = matSxC$Sample
        matSxC$Sample = NULL
        matSxC = as.matrix(matSxC)

        # Should be sorted but just to be sure
        matSxC = matSxC[, order(thisModel$Mean) ]
        colnames(matSxC) = paste0( thisFeature, 1:ncol(matSxC) )

        return(matSxC)

    } )

    # Return data
    allMats = do.call(cbind, lMats)
    lMats = list(sampleByComponent = allMats, model=allModels)
    return(lMats)
}
