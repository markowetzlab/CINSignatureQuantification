#' estimateThresholds
#'
#' @param object An object of class SigQuant or PLACEHOLDERCLASS containing computed features, signature definitions, siganture activities, and feature mixture models.
#' @param iters Number of noise simulation iterations to perform to estimate thresholds.
#' @param method Feature and mixture model method to use for computing simulated features and signature activities
#' @param parallel Use doFuture parallelisation to run iterations set by future::plan()
#'
#' @return data.frame containing series of signature-specific thresholds (Default used = Thresh_ZeroGMM_0.05)
#' @import mclust
#' @export estimateThresholds
#'
estimateThresholds <- function(object=NULL,iters=1000,SDPROP = 20,FINALWIDTH = 0.1,RANGECNAS = 0.1,
                               method="drewsV2",minmaxNorm=FALSE,orderOri=FALSE,parallel=FALSE){
    if(is.null(object)){
        stop("no data")
    }
    if(iters < 2){
        stop("require at least 2 iterations")
    }

    feats <- object@featData
    sigDefs <- object@backup.signatures
    sigActs <- object@activities$normAct1
    mixtures <- object@featFitting$model

    if(parallel) {
        if (!requireNamespace("doFuture", quietly = TRUE)) {
            stop(
                "Package \"doFuture\" must be installed to use multiple threads/cores.",
                call. = FALSE
            )
        }

        message(paste0("running doFuture foreach for ",iters," iterations"))
        # Multi-thread usage
        `%dofuture%` <- doFuture::`%dofuture%`
        # provide progressr call if wanted
        p <- progressr::progressor(steps = iters)
        lActsList <- foreach::foreach(i=1:iters,.options.future = list(seed = TRUE,packages = c("CINSignatureQuantification"))) %dofuture% {
            #message(paste0("Noise sim iteration ",i," of ",iters))
            lFeatures <- simulateFeatureNoise(feats=feats,SDPROP=SDPROP,FINALWIDTH=FINALWIDTH,RANGECNAS=RANGECNAS)
            switch(method,
                   drews={
                       lSxC <- calculateSampleByComponentMatrixDrews(brECNF = lFeatures,
                                                                     UNINFPRIOR = TRUE)$sampleByComponent
                       acts <- calculateActivityDrews(object = lSxC)
                       lActs = t(apply(acts, 2, function(x) x/sum(x)))
                   },
                   drewsV2={
                       lSxC <- calculateSampleByComponentMatrixDrewsV2(brECNF = lFeatures,
                                                                       UNINFPRIOR = TRUE)$sampleByComponent
                       ## to be changed once V2 is updated with new definitions
                       acts <- calculateActivityDrewsV2(object = lSxC)
                       lActs = t(apply(acts, 2, function(x) x/sum(x)))
                   })
            # progressr iter call
            p(sprintf("i=%g", i))
            list(lActs)
        }
    } else {
        lActsList <- list()
        for(i in 1:iters){
            message(paste0("running sequential iterations for  ",i," of ",iters," iterations"))
            lFeatures <- simulateFeatureNoise(feats=feats,SDPROP=SDPROP,FINALWIDTH=FINALWIDTH,RANGECNAS=RANGECNAS)

            switch(method,
                    drews={
                        lSxC <- calculateSampleByComponentMatrixDrews(brECNF = lFeatures,
                                                                      UNINFPRIOR = TRUE)$sampleByComponent
                        acts <- calculateActivityDrews(object = lSxC)
                        lActs = t(apply(acts, 2, function(x) x/sum(x)))
                    },
                    drewsV2={
                        lSxC <- calculateSampleByComponentMatrixDrewsV2(brECNF = lFeatures,
                                                                        UNINFPRIOR = TRUE)$sampleByComponent
                        ## to be changed once V2 is updated with new definitions
                        acts <- calculateActivityDrewsV2(object = lSxC)
                        lActs = t(apply(acts, 2, function(x) x/sum(x)))
            })
            lActsList <- append(lActsList,list(lActs))
        }
    }
    message(paste0("Computing thresholds from ",iters, " simulated Signature activities"))
    thresholdTab <- calculateThresholds(lSignatures = lActsList,originalActivities = sigActs,
                        minmaxNorm = minmaxNorm,orderOri = orderOri)
    thresholdTab <- thresholdTab[rownames(sigDefs),]
    return(thresholdTab)
}

simulateFeatureNoise <- function(feats=NULL,SDPROP = 20,FINALWIDTH = 0.1,RANGECNAS = 0.1){

    lFeatures = sapply(names(feats), function(thisFeat) {

        dfFeat = feats[[thisFeat]]

        if(thisFeat %in% c("segsize", "changepoint", "copynumber")) {
            # Pass through option
            if(FINALWIDTH != 0) {
                dfFeat = addGaussianNoise(dfFeat, sdProp = SDPROP, finalWidth = FINALWIDTH)
            }
        } else {
            # Pass through option
            if(RANGECNAS != 0) {
                dfFeat = addSampleNoise(dfFeat, RANGECNAS = RANGECNAS)
            }
        }
        return(dfFeat)
    }, simplify = FALSE, USE.NAMES = TRUE)
}

addSampleNoise <- function(dfFeat, RANGECNAS = 0.1) {

    # bpchrarm has a different column name
    colnames(dfFeat) = c("ID", "value")

    # Split into sample-specific lists for easier looping
    dfFeat$ID = factor(dfFeat$ID, levels = unique(dfFeat$ID))
    lFeat = split(dfFeat, dfFeat$ID)
    lSim = lapply(lFeat, function(thisSample) {

        numCNAs = sum(thisSample$value)
        tabSample = table(thisSample$value)

        # Shortcut for osCN feature where vectors of only zeros is possible
        if(numCNAs == 0) { return(thisSample) }

        # Extreme case if only one non-zero entry is present (apparently sample function fails)
        # Implement rescue to original values if only one value is present (independent of how many)
        if(length(tabSample) == 1) { return(thisSample) }

        # Draw new values with probabilities of existing sample
        vSim = sample(x = as.numeric(names(tabSample)), size = sum(tabSample),
                      prob = tabSample/sum(tabSample), replace = TRUE)
        currentDraw = sum(vSim)

        # Check if this draw is above the user-defined limit. If so, enter while loop and draw until resolved
        while( abs(1-(currentDraw/numCNAs)) > RANGECNAS) {
            vSim = sample(x = as.numeric(names(tabSample)), size = sum(tabSample),
                          prob = tabSample/sum(tabSample), replace = TRUE)
            currentDraw = sum(vSim)
        }

        thisSample$value = vSim
        return(thisSample)
    })

    # Merge list and prepare for return
    dfSim = do.call(rbind, lSim)
    rownames(dfSim) = NULL
    dfSim$ID = as.character(dfSim$ID)

    return(dfSim)
}

addGaussianNoise <- function(dfFeat, sdProp = 20, finalWidth = NA) {

    ## This function adds noise to a column called "value".
    ## The level of noise is determined by the variable "sdProp" - sd = value / sdProp
    ## If wished the width of noise can be limited by "finalWidth" to a determined maximum width around
    ## the original value (to avoid large outliers).

    dfFeat$value2 = rnorm(n = length(dfFeat$value), mean = dfFeat$value, sd = dfFeat$value/sdProp)

    # Remove outliers if wished
    if(! is.na(finalWidth)) {

        dfFeat$diff = dfFeat$value / dfFeat$value2
        dfFeat$diff[is.na(dfFeat$diff)] <- 1
        dfFeat$redraw = abs(1-dfFeat$diff) > finalWidth
        while(sum(dfFeat$redraw) > 0) {

            dfFeat$value2[dfFeat$redraw] = stats::rnorm(n = length(dfFeat$value[dfFeat$redraw]),
                                                 mean = dfFeat$value[dfFeat$redraw],
                                                 sd = dfFeat$value[dfFeat$redraw]/sdProp)
            dfFeat$diff = dfFeat$value / dfFeat$value2
            dfFeat$diff[is.na(dfFeat$diff)] <- 1
            dfFeat$redraw = abs(1-dfFeat$diff) > finalWidth
        }

    }

    # Clean up and return
    dfFeat$diff = NULL
    dfFeat$redraw = NULL
    dfFeat$value = dfFeat$value2
    dfFeat$value2 = NULL

    return(dfFeat)
}

calculateThresholds <- function(lSignatures=NULL,originalActivities=NULL,minmaxNorm=FALSE,orderOri=FALSE) {

    # Bring list of signature activities into one data table
    dtSigs <- data.table::data.table(formatSignatureList(lSignatures),stringsAsFactors = TRUE)
    siglevels <- stringr::str_sort(levels(factor(dtSigs$Var2)),numeric = T)
    dtSigs$Var2 <- factor(dtSigs$Var2,levels = siglevels)

    dtOri <- data.table(formatSignatureList(originalActivities))
    dtOri$Var2 <- factor(dtOri$Var2,levels = siglevels)

    # Plot boxplot for each sample and signature
    allSigs = siglevels
    lthresh = lapply(allSigs, function(thisSig) {

        dtCS2 = dtSigs[dtSigs$Var2 == thisSig,]

        # prepare original values
        dtOriCS2 = dtOri[ dtOri$Var2 == thisSig,]

        # Sort samples by original signature value or by median simulation value
        if(orderOri) {
            newOrder = as.character(dtOriCS2$Var1[ order(dtOriCS2$value, decreasing = TRUE)])
        } else {
            dtOrder = stats::aggregate(value ~ Var1, dtCS2, stats::median)
            newOrder = as.character(dtOrder$Var1[ order(dtOrder$value, decreasing = TRUE)])
        }
        # Reorder factors
        dtOriCS2$Var1 = factor(as.character(dtOriCS2$Var1), levels = newOrder)
        dtCS2$Var1 = factor(as.character(dtCS2$Var1), levels = newOrder)

        # Reorder the names themselves (needed for identifying the correct threshold later)
        ## fixed and replaced match(newOrder, dtOri$Var1) whihc returned NA.table
        dtOriCS2 = dtOriCS2[ match(newOrder, dtOriCS2$Var1), ]

        ## Four methods for identifying thresholds
        # Method 1: First time an IQR crosses zero
        dfQ25 = aggregate(value ~ Var1, dtCS2, q25)
        # Order of samples should match the sorting from above aka the factor levels
        if(!identical(as.character(dfQ25$Var1),levels(dfQ25$Var1))){
            dfQ25 = dfQ25[ match(newOrder, dfQ25$Var1), ]
        }
        # Identify first index, corresponding sample and original value when simulation reached zero
        #firstZeroIndex = rle(dfQ25$value <= 0)$lengths[1]+1 # not sure why rle is used
        firstZeroIndex <- min(which(dfQ25$value == 0))
        thresh1 = dtOriCS2$value[firstZeroIndex]

        # Method 2: 95% quantile of all true zero samples
        allZeroSamples = as.character(dtOriCS2$Var1[ dtOriCS2$value == 0 ])
        dfQ75 = aggregate(value ~ Var1, dtCS2[ dtCS2$Var1 %in% allZeroSamples, ], q75)
        thresh2 = q95(dfQ75$value)
        # Method 3 & 4: Use values from true zero samples to fit a Gaussian mixture model and identify
        # first samples that don't fit into the Gaussian.
        # Threshold 3: <5%
        # Threshold 4: <1%
        dat = dtCS2$value[ dtCS2$Var1 %in% as.character(dtOriCS2$Var1[ dtOriCS2$value == 0 ]) ]
        gModel = mclust::Mclust(dat, modelNames = "V", G = 1, verbose = FALSE)
        modelMean = gModel$parameters$mean
        modelVar = gModel$parameters$variance$sigmasq
        thresh3 = stats::qnorm(0.95, mean = modelMean, sd = sqrt(modelVar), lower.tail = TRUE)
        thresh4 = stats::qnorm(0.99, mean = modelMean, sd = sqrt(modelVar), lower.tail = TRUE)

        return(c(thresh1, thresh2, thresh3, thresh4))
    })

    names(lthresh) = allSigs
    lthresh <- do.call(rbind,lapply(lthresh,FUN = function(x) x))
    colnames(lthresh) = c("Thresh_Q25ZeroHit", "Thresh_Q95ofZeroQ75",
                            "Thresh_ZeroGMM_0.05", "Thresh_ZeroGMM_0.01")
    return(lthresh)
}

formatSignatureList <- function(list){

    if(!is.null(dim(list)) & is.matrix(list)){
        y <- as.data.frame(list)
        y$Var1 <- rownames(y)
        y<- tidyr::pivot_longer(y,cols = 1:ncol(y)-1)
        colnames(y) <- c("Var1","Var2","value")
        formattedList <- data.frame(y)
    } else {
        formattedList <- do.call(rbind,lapply(1:length(list),FUN = function(x){
            y <- as.data.frame(list[[x]])
            y$Var1 <- rownames(y)
            y<- tidyr::pivot_longer(y,cols = 1:ncol(y)-1)
            y$iter <- x
            colnames(y) <- c("Var1","Var2","value","iter")
            return(data.frame(y))
        }))
    }
    return(formattedList)
}

q25 <- function(x){
    stats::quantile(x, probs = 0.25)
}
q75 <- function(x){
    stats::quantile(x, probs = 0.75)
}
q95 <- function(x){
    stats::quantile(x, probs = 0.95)
}
