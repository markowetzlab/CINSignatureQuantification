# estimateActivityThresholds

library(data.table)
library(YAPSA)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(Cairo)
library(mclust)

OUTPATH=file.path("data/MC_simulation_overfitting")
SAVEINTERMEDIATEFILES=TRUE
# Misc id
NAMEID="_fullTCGA_1000sims_10pGaussian_10pSamplePoisson"
## Feats
# ORIECNF=file.path(BASE, "out/2_TCGA_PCAWG_ECNF.rds")
# lOriECNF = readRDS(ORIECNF)
lOriECNF <- tcga_6335_cnquant_newweights@featData
# Simulations
NUMSIMS=1e3
## mixtures
INPUTMODELS=file.path(BASE, "CINsignatures/Mixmodels_merged_components.rds")
allModels = readRDS(INPUTMODELS)
UNINFPRIOR=TRUE
# SigDefs
# SIGNATUREFILE=file.path(BASE, "CINsignatures/Signature_Compendium_v5_Cosine-0.74_Signatures_NAMESAPRIL21.rds")
# W = readRDS(SIGNATUREFILE)
W = matAll[rownames(matAll) %in% colnames(compOldNewCosineFiltAll_074_sub),]
# SigActs
# ACTIVITIES=file.path(BASE, "CINsignatures/Signature_Compendium_v5_Cosine-0.74_Activities_NAMESAPRIL21.rds")
# dtOri = data.table(melt(readRDS(ACTIVITIES)))
dtOri <- data.table(melt(lSigs))
## Step 1: Simulate feature distributions with noise
addSampleNoise = function(dfFeat, RANGECNAS = 0.1) {

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

# Segment size, changepoint => Gaussian noise
addGaussianNoise = function(dfFeat, sdProp = 20, finalWidth = NA) {

    ## This function adds noise to a column called "value".
    ## The level of noise is determined by the variable "sdProp" - sd = value / sdProp
    ## If wished the width of noise can be limited by "finalWidth" to a determined maximum width around
    ## the original value (to avoid large outliers).

    dfFeat$value2 = rnorm(n = length(dfFeat$value), mean = dfFeat$value, sd = dfFeat$value/sdProp)

    # Remove outliers if wished
    if(! is.na(finalWidth)) {

        dfFeat$diff = dfFeat$value / dfFeat$value2
        dfFeat$redraw = abs(1-dfFeat$diff) > finalWidth
        while(sum(dfFeat$redraw) > 0) {

            dfFeat$value2[dfFeat$redraw] = rnorm(n = length(dfFeat$value[dfFeat$redraw]),
                                                 mean = dfFeat$value[dfFeat$redraw],
                                                 sd = dfFeat$value[dfFeat$redraw]/sdProp)
            dfFeat$diff = dfFeat$value / dfFeat$value2
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

addNoiseToFeatures = function(lOriECNF, allFeatures = c("segsize", "changepoint", "bp10MB", "osCN", "bpchrarm"),
                              SDPROP = 20, FINALWIDTH = 0.1, RANGECNAS = 0.1, NUMSIMS = 1) {

    lSim = lapply(1:NUMSIMS, function(thisN) {

        if(thisN %% 10 == 0) {
            print(thisN)
        }
        lFeatures = sapply(allFeatures, function(thisFeat) {

            dfFeat = lOriECNF[[thisFeat]]

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

    } )

    return(lSim)
}

lSimulation = addNoiseToFeatures(lOriECNF, allFeatures = c("changepoint", "segsize",
                                                           "bpchrarm", "osCN", "bp10MB"),
                                 SDPROP = 20, FINALWIDTH = 0.1, RANGECNAS = 0.1, NUMSIMS = 10)


## Step 2: Derive SxC matrices
print("Start deriving SxC matrices data...")
# lMatrices = deriveSxCMatrices(lSimulation, allModels = allModels,
#                               allFeatures = names(allModels), UNINFPRIOR = TRUE)
lMatrices <- lapply(X = lSimulation,FUN = function(x)
    CINSignatureQuantification:::calculateSampleByComponentMatrixDrewsV2(brECNF = x,UNINFPRIOR = TRUE)$sampleByComponent)

# if(SAVEINTERMEDIATEFILES) {
#     saveRDS(lMatrices, file.path(OUTPATH, paste0("1_SxC_Matrices", NAMEID, ".rds")))
# }

## Step 3: Calculate signature activities
print("Start calculating signature activities...")

# lSignatures = calculateSignatureActivities(lMatrices, W = W)
lSignatures <- lapply(lMatrices,FUN = function(x){
    y <- CINSignatureQuantification:::calculateActivityDrewsV2(x, SIGS = W)
    H = t( apply(y, 2, function(x) x/sum(x)) )
    return(H)
})


estimateThresholds <- function(feats=NULL,sigDefs=NULL,sigActs=NULL,mixtures=NULL,
                               iters=1000,SDPROP = 20,FINALWIDTH = 0.1, RANGECNAS = 0.1,UNINFPRIOR = TRUE,method=NULL){
    lActsList <- list()
    for(i in 1:iters){

        lFeatures <- simulateFeatureNoise(feats=feats,SDPROP=SDPROP,FINALWIDTH=FINALWIDTH,RANGECNAS=RANGECNAS)

        switch(method,
                drews={
                    lSxC <- calculateSampleByComponentMatrixDrews(brECNF = lFeatures,
                                                                  UNINFPRIOR = TRUE)$sampleByComponent
                    acts <- CINSignatureQuantification:::calculateActivityDrews(object = lSxC)
                    lActs = t(apply(acts, 2, function(x) x/sum(x)))
                },
                drewsV2={
                    lSxC <- calculateSampleByComponentMatrixDrewsV2(brECNF = lFeatures,
                                                                    UNINFPRIOR = TRUE)$sampleByComponent
                    ## to be changed once V2 is updated with new definitions
                    acts <- CINSignatureQuantification:::calculateActivityDrewsV2(object = lSxC,SIGS = sigDefs)
                    lActs = t(apply(acts, 2, function(x) x/sum(x)))
        })
        lActsList <- append(lActsList,lActs)
    }
    return(lActsList)
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

addSampleNoise = function(dfFeat, RANGECNAS = 0.1) {

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
# Segment size, changepoint => Gaussian noise
addGaussianNoise = function(dfFeat, sdProp = 20, finalWidth = NA) {

    ## This function adds noise to a column called "value".
    ## The level of noise is determined by the variable "sdProp" - sd = value / sdProp
    ## If wished the width of noise can be limited by "finalWidth" to a determined maximum width around
    ## the original value (to avoid large outliers).

    dfFeat$value2 = rnorm(n = length(dfFeat$value), mean = dfFeat$value, sd = dfFeat$value/sdProp)

    # Remove outliers if wished
    if(! is.na(finalWidth)) {

        dfFeat$diff = dfFeat$value / dfFeat$value2
        dfFeat$redraw = abs(1-dfFeat$diff) > finalWidth
        while(sum(dfFeat$redraw) > 0) {

            dfFeat$value2[dfFeat$redraw] = rnorm(n = length(dfFeat$value[dfFeat$redraw]),
                                                 mean = dfFeat$value[dfFeat$redraw],
                                                 sd = dfFeat$value[dfFeat$redraw]/sdProp)
            dfFeat$diff = dfFeat$value / dfFeat$value2
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

q25 = function(x) { quantile(x, probs = 0.25)}
q75 = function(x) { quantile(x, probs = 0.75)}
q95 = function(x) { quantile(x, probs = 0.95)}

###
tcga_6335_cnquant_newweights <- readRDS("data/tcga_6335_cnquant_newweights.rds")
lOriECNF <- getFeatures(tcga_6335_cnquant_newweights)
W <- sigDef
lSigs <- lSigs
INPUTMODELS <- tcga_6335_cnquant_newweights@featFitting$model

thresholdActs <- estimateThresholds(feats = lOriECNF,sigDefs = W,sigActs = lSigs,mixtures = INPUTMODELS,iters = 2)


# if(SAVEINTERMEDIATEFILES) {
#     saveRDS(lSignatures, file.path(OUTPATH, paste0("2_Activities", NAMEID, ".rds")))
# }

## Plot boxplot for each signature and sample
plotAllSigs = function(lSignatures, dtOri, minmaxNorm = FALSE, orderOri = FALSE, DOTSIZE = 0.25,
                       TITLE = "478 TCGA/PCAWG samples", FILEOUT = NULL) {

    # Check if output path exists
    if(is.null(FILEOUT)) {
        warning("No output path for threshold text file supplied. Defaulting to './'!")
        #FILEOUT = "./Thresholds.txt"
    }

    # Bring list of signature activities into one data table
    lMelt = list()
    allSigs = length(lSignatures)
    for(i in 1:allSigs) {
        dtSignatures = data.table(melt(lSignatures[[i]] ))
        dtSignatures$iter = i
        lMelt[[length(lMelt)+1]] = dtSignatures
    }

    dtSigs = rbindlist(lMelt)

    # Plot boxplot for each sample and signature
    allSigs = levels(dtSigs$Var2)
    lPlots = lapply(allSigs, function(thisSig) {

        dtCS2 = dtSigs[ dtSigs$Var2 == thisSig, ]

        # prepare original values
        dtOriCS2 = dtOri[ dtOri$Var2 == thisSig, ]

        # Sort samples by original signature value or by median simulation value
        if(orderOri) {
            newOrder = as.character(dtOriCS2$Var1[ order(dtOriCS2$value, decreasing = TRUE) ])
        } else {
            dtOrder = aggregate(value ~ Var1, dtCS2, median)
            newOrder = as.character(dtOrder$Var1[ order(dtOrder$value, decreasing = TRUE) ])
        }
        # Reorder factors
        dtOriCS2$Var1 = factor(as.character(dtOriCS2$Var1), levels = newOrder)
        dtCS2$Var1 = factor(as.character(dtCS2$Var1), levels = newOrder)

        # Reorder the names themselves (needed for identifying the correct threshold later)
        dtOriCS2 = dtOriCS2[ match(newOrder, dtOri$Var1), ]

        ## Four methods for identifying thresholds
        # Method 1: First time an IQR crosses zero
        dfQ25 = aggregate(value ~ Var1, dtCS2, q25)
        # Order of samples should match the sorting from above aka the factor levels
        if(! identical(as.character(dfQ25$Var1), levels(dfQ25$Var1))) { dfQ25 = dfQ25[ match(newOrder, dfQ25$Var1), ] }
        # Identify first index, corresponding sample and original value when simulation reached zero
        firstZeroIndex = rle(dfQ25$value <= 0)$lengths[1]+1
        thresh1 = dtOriCS2$value[ firstZeroIndex ]


        # Method 2: 95% quantile of all true zero samples
        allZeroSamples = as.character(dtOriCS2$Var1[ dtOriCS2$value == 0 ])

        dfQ75 = aggregate(value ~ Var1, dtCS2[ dtCS2$Var1 %in% allZeroSamples, ], q75)
        thresh2 = q95(dfQ75$value)


        # Method 3 & 4: Use values from true zero samples to fit a Gaussian mixture model and identify
        # first samples that don't fit into the Gaussian.
        # Threshold 3: <5%
        # Threshold 4: <1%
        dat = dtCS2$value[ dtCS2$Var1 %in% as.character(dtOriCS2$Var1[ dtOriCS2$value == 0 ]) ]
        gModel = Mclust(dat, modelNames = "V", G = 1, verbose = FALSE)
        modelMean = gModel$parameters$mean
        modelVar = gModel$parameters$variance$sigmasq
        thresh3 = qnorm(0.95, mean = modelMean, sd = sqrt(modelVar), lower.tail = TRUE)
        thresh4 = qnorm(0.99, mean = modelMean, sd = sqrt(modelVar), lower.tail = TRUE)

        # Plot IQR
        p1 = ggplot(dtCS2, aes(x = Var1, y = value)) + geom_boxplot(outlier.shape = NA, coef = 0) +
            #geom_point(data = dtOriCS2, aes(x = Var1, y = value), colour = "red", size = DOTSIZE) +
            #geom_point(aes(y = ori), colour = "red", size = DOTSIZE) +
            geom_hline(yintercept = thresh1, size = 1, linetype = "dotted", colour = "red") +
            geom_vline(xintercept = firstZeroIndex, size = 1, linetype = "dotted", colour = "red") +
            geom_hline(yintercept = thresh2, size = 1, linetype = "dashed", colour = "blue") +
            geom_hline(yintercept = thresh3, size = 1, linetype = "dotdash", colour = "green") +
            geom_hline(yintercept = thresh4, size = 1, linetype = "dotdash", colour = "darkgreen") +
            #annotate("text", x = 400, y = thresh1*1.5, size = 6, label = as.character(signif(thresh1, digits = 2))) +
            #annotate("text", x = 400, y = 0, size = 6, label = as.character(signif(thresh2, digits = 2))) +
            #annotate("text", x = 6000, y = 0, size = 6, label = as.character(signif(thresh3, digits = 2))) +
            #annotate("text", x = 6000, y = thresh4*1.5, size = 6, label = as.character(signif(thresh4, digits = 2))) +
            theme_bw() +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank()) +
            ggtitle(TITLE) + ylab(paste(thisSig)) +
            scale_x_discrete(breaks = NULL)

        return(list(p1, c(thresh1, thresh2, thresh3, thresh4)))
    })

    names(lPlots) = allSigs

    # Split plots and thresholds
    lOutPlot = list()
    lOutThresh = list()
    for(thisSig in allSigs) {
        lOutPlot[[length(lOutPlot)+1]] = lPlots[[thisSig]][[1]]
        lOutThresh[[length(lOutThresh)+1]] = c(thisSig, signif(lPlots[[thisSig]][[2]], digits = 4))
    }

    # Combine thresholds and save output file
    dtOutThresh = data.table(do.call(rbind, lOutThresh))
    colnames(dtOutThresh) = c("Sig", "Thresh_Q25ZeroHit", "Thresh_Q95ofZeroQ75",
                              "Thresh_ZeroGMM_0.05", "Thresh_ZeroGMM_0.01")
    #write.table(dtOutThresh, FILEOUT, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

    return(dtOutThresh)
}

library(mclust)
allSigs <- plotAllSigs(lSignatures, dtOri, orderOri = TRUE, DOTSIZE = 0.25, TITLE = "1000 runs full TCGA")
