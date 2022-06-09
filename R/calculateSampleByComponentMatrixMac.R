calculateSampleByComponentMatrixMac = function(brECNF, UNINFPRIOR = FALSE) {
    # This version of the SampleByComponent calculation uses the same
    # methodology as Drews 2022 to calculate Macintyre 2018 SxC matrix without
    # the need for flexmix. This includes the model weighting (UNIFPRIOR=FALSE)
    # which is in contrast to the Drews method (UNIFPRIOR=TRUE). The
    # corresponding SxC matrix is not identical to the original Macintyre 2018
    # method with a mean and median difference in sum of posteriors of -3.33e-13
    # & -2.12e-14 across 478 TCGA gold standard samples.

    # Load mix models
    #allModels = get(load("data/Macintyre2018_OV_Mixture_Models.rda"))
    allModels = get(data("Macintyre2018_OV_Mixture_Models",envir = environment()))
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

## OLD method for generating SxC using flexmix
# calculateSampleByComponentMatrixMacOLD <- function(CN_features, all_components=NULL, cores = 1, rowIter = 1000, subcores = 2, method=method)
# {
#     if(is.null(all_components))
#     {
#         # This file is large and should be reimplemented in some form - Specifically only retain data flexmix::posterior requires for SoP?
#         all_components<-get(load("inst/mac_component_parameters.rda"))
#     }
#     # if(cores > 1){
#     #     require(foreach)
#     #     feats = c( "segsize", "bp10MB", "osCN", "changepoint", "copynumber", "bpchrarm" )
#     #     doMC::registerDoMC(cores)
#     #     full_mat = foreach(feat=feats, .combine=cbind) %dopar% {
#     #         calculateSumOfPosteriors_MAC(CN_features[[feat]],all_components[[feat]],
#     #                                  feat, rowIter = rowIter, cores = subcores)
#     #     }
#     # } else {
#     full_mat<-cbind(
#         calculateSumOfPosteriors_MAC(CN_features[["segsize"]],all_components[["segsize"]],"segsize"),
#         calculateSumOfPosteriors_MAC(CN_features[["bp10MB"]],all_components[["bp10MB"]],"bp10MB"),
#         calculateSumOfPosteriors_MAC(CN_features[["osCN"]],all_components[["osCN"]],"osCN"),
#         calculateSumOfPosteriors_MAC(CN_features[["changepoint"]],all_components[["changepoint"]],"changepoint"),
#         calculateSumOfPosteriors_MAC(CN_features[["copynumber"]],all_components[["copynumber"]],"copynumber"),
#         calculateSumOfPosteriors_MAC(CN_features[["bpchrarm"]],all_components[["bpchrarm"]],"bpchrarm"))
#     #}
#     rownames(full_mat)<-unique(CN_features[["segsize"]][,1])
#     full_mat[is.na(full_mat)]<-0
#     list(method=method,sampleByComponent=full_mat,model=all_components)
# }
#
# calculateSumOfPosteriors_MAC<-function(CN_feature,components,name, rowIter = 1000, cores = 1)
# {
#     # if(cores > 1){
#     #     require(foreach)
#     #     require(doMC)
#     #     len = dim(CN_feature)[1]
#     #     iters = floor( len / rowIter )
#     #     lastiter = iters[length(iters)]
#     #     registerDoMC(cores)
#     #     curr_posterior = foreach( i=0:iters, .combine=rbind) %dopar% {
#     #         start = i*rowIter+1
#     #         if(i != lastiter) { end = (i+1)*rowIter } else { end = len }
#     #         flexmix::posterior(components,data.frame(dat=as.numeric(CN_feature[start:end,2])))
#     #     }
#     # } else {
#     curr_posterior<-flexmix::posterior(components,data.frame(dat=as.numeric(CN_feature[,2])))
#     #}
#     mat<-cbind(CN_feature,curr_posterior)
#     posterior_sum<-c()
#     ## foreach and parallelising doesn't make the following code faster.
#     for(i in unique(mat$ID))
#     {
#         posterior_sum<-rbind(posterior_sum,colSums(mat[mat$ID==i,c(-1,-2)]))
#     }
#     params<-flexmix::parameters(components)
#     if(!is.null(nrow(params)))
#     {
#         posterior_sum<-posterior_sum[,order(params[1,])]
#     }
#     else
#     {
#         posterior_sum<-posterior_sum[,order(params)]
#     }
#     colnames(posterior_sum)<-paste0(name,1:ncol(posterior_sum))
#     rownames(posterior_sum)<-rownames(unique(mat$ID))
#     posterior_sum
# }
