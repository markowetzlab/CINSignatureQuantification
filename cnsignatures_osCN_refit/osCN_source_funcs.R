# osCN feature testing functions
## condensed functions used in testing feature fix
#### Functions
## original pre-package functions
getOscillationDrews = function(abs_profiles, chrlen) {

    # Prepare looping
    out = c()
    samps = names(abs_profiles)
    # Loop over samples
    for(i in samps) {

        # Retrieve segments
        segTab = abs_profiles[[i]]
        colnames(segTab)[4] = "segVal"

        # Loop over chromosomes to identify oscillation
        chrs = unique(segTab$chromosome)
        oscCounts = c()
        for(c in chrs) {

            currseg = as.numeric(segTab$segVal[segTab$chromosome == c])
            currseg = round(as.numeric(currseg))

            # Only take chains into consideration with a length of more than 3 elements
            if(length(currseg)>3) {
                prevval = currseg[1]
                count = 0
                for(j in 3:length(currseg)) {
                    if(currseg[j] == prevval & currseg[j] != currseg[j-1]) {
                        count = count+1
                        # # suggested fix for missing end of chromosome chain counts
                        # # https://github.com/markowetzlab/CINSignatureQuantification/issues/22
                        # if (j == length(currseg)) {
                        #     oscCounts = c(oscCounts, count)
                        #     count = 0
                        # }
                    } else {
                        oscCounts = c(oscCounts,count)
                        count = 0
                    }
                    prevval = currseg[j-1]
                }
            }
        }
        # Make sure it's really numeric
        out = rbind(out, cbind(ID = rep(i,length(oscCounts)),
                               value = as.numeric(oscCounts)))
        if(length(oscCounts) == 0) {
            out = rbind(out,cbind(ID = i, value = 0))
        }
    }

    # Prepare return
    rownames(out) = NULL
    return(data.frame(out,stringsAsFactors = FALSE))
}

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
        ## Added fix for circumstances where zero osCN features are recorded
        ## Should also be applicable to other features - CINSignaturesQuantification issue #17
        if(is.null(dim(postDatUnscaled))){
            postDatScaled = data.frame(t(postDatUnscaled / sum(postDatUnscaled)))
        } else {
            postDatScaled = data.frame(postDatUnscaled / rowSums(postDatUnscaled))
        }
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

## pre-package V2 functions
getOscillationDrewsV2 = function(abs_profiles, chrlen) {

    # Prepare looping
    out = c()
    samps = names(abs_profiles)
    # Loop over samples
    for(i in samps) {

        # Retrieve segments
        segTab = abs_profiles[[i]]
        colnames(segTab)[4] = "segVal"

        # Loop over chromosomes to identify oscillation
        chrs = unique(segTab$chromosome)
        oscCounts = c()
        for(c in chrs) {

            currseg = as.numeric(segTab$segVal[segTab$chromosome == c])
            currseg = round(as.numeric(currseg))

            # Only take chains into consideration with a length of more than 3 elements
            if(length(currseg)>3) {
                prevval = currseg[1]
                count = 0
                for(j in 3:length(currseg)) {
                    if(currseg[j] == prevval & currseg[j] != currseg[j-1]) {
                        count = count+1
                        # suggested fix for missing end of chromosome chain counts
                        # https://github.com/markowetzlab/CINSignatureQuantification/issues/22
                        if (j == length(currseg)) {
                            oscCounts = c(oscCounts, count)
                            count = 0
                        }
                    } else {
                        oscCounts = c(oscCounts,count)
                        count = 0
                    }
                    prevval = currseg[j-1]
                }
            }
        }
        # Make sure it's really numeric
        out = rbind(out, cbind(ID = rep(i,length(oscCounts)),
                               value = as.numeric(oscCounts)))
        if(length(oscCounts) == 0) {
            out = rbind(out,cbind(ID = i, value = 0))
        }
    }

    # Prepare return
    rownames(out) = NULL
    return(data.frame(out,stringsAsFactors = FALSE))
}

## pre-package V2 functions with labelling of chain location
getOscillationDrewsV2WITHLABEL = function(abs_profiles, chrlen) {

    # Prepare looping
    out = c()
    samps = names(abs_profiles)
    # Loop over samples
    for(i in samps) {

        # Retrieve segments
        segTab = abs_profiles[[i]]
        colnames(segTab)[4] = "segVal"

        # Loop over chromosomes to identify oscillation
        chrs = unique(segTab$chromosome)
        oscCounts = c()
        osLabel = c()
        for(c in chrs) {

            currseg = as.numeric(segTab$segVal[segTab$chromosome == c])
            currseg = round(as.numeric(currseg))

            # Only take chains into consideration with a length of more than 3 elements
            if(length(currseg)>3) {
                prevval = currseg[1]
                count = 0
                for(j in 3:length(currseg)) {
                    if(currseg[j] == prevval & currseg[j] != currseg[j-1]) {
                        count = count+1
                        # suggested fix for missing end of chromosome chain counts
                        # https://github.com/markowetzlab/CINSignatureQuantification/issues/22
                        if (j == length(currseg)) {
                            oscCounts = c(oscCounts, count)
                            osLabel = c(osLabel, TRUE)
                            count = 0
                        }
                    } else {
                        oscCounts = c(oscCounts,count)
                        osLabel = c(osLabel, FALSE)
                        count = 0
                    }
                    prevval = currseg[j-1]
                }
            }
        }
        # Make sure it's really numeric
        out = rbind(out, cbind(ID = rep(i,length(oscCounts)),
                               value = as.numeric(oscCounts),
                                label = as.logical(osLabel)))
        if(length(oscCounts) == 0) {
            out = rbind(out,cbind(ID = i, value = 0,label = FALSE))
        }
    }

    # Prepare return
    rownames(out) = NULL
    return(data.frame(out,stringsAsFactors = FALSE))
}


calculateSampleByComponentMatrixDrewsV2 = function(brECNF, UNINFPRIOR = TRUE,weights=Drews2022_v2_TCGA_Mixture_Models) {

    # Load mix models
    #allModels = get(load("data/Drews2022_TCGA_Mixture_Models.rda"))
    allModels = Drews2022_v2_TCGA_Mixture_Models
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
        ## Added fix for circumstances where zero osCN features are recorded
        ## Should also be applicable to other features - CINSignaturesQuantification issue #17
        if(is.null(dim(postDatUnscaled))){
            postDatScaled = data.frame(t(postDatUnscaled / sum(postDatUnscaled)))
        } else {
            postDatScaled = data.frame(postDatUnscaled / rowSums(postDatUnscaled))
        }
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

calculateActivityDrewsV2 = function(object,cancer.subset=NULL,SIGS=NULL) {

    # Extract relevant information from object
    # V = object@featFitting$sampleByComponent
    # nSamp = nrow(object@featFitting$sampleByComponent)
    # nFeat = ncol(object@featFitting$sampleByComponent)
    V = object
    nSamp = nrow(object)
    nFeat = ncol(object)


    # Load signatures
    #W = get(load("data/Drews2022_TCGA_Signatures.rda"))
    W = SIGS
    if(!is.null(cancer.subset)){
        subset = getCancerSpecificSignatures(cancer.subset)
        W = W[rownames(W) %in% subset,]
    }
    # Sanity check mutational catalogue (not really necessary)
    if(nSamp > nFeat) {
        # Case 1: More samples than features
        if(nrow(V) > ncol(V)) { V = t(V) }
    } else if(nSamp < nFeat) {
        # Case 2: Fewer samples than features
        if(nrow(V) < ncol(V)) { V = t(V) }
    } else {
        # Case 3: Edge case where there are as many samples as features
        if(sum(grepl("segsize", colnames(V))) > 0) { V = t(V) }
    }


    # Sanity check signature matrix
    if(nrow(W) < ncol(W)) W = t(W)
    # Check order of components and fix if necessary
    if(! identical(rownames(W), rownames(V)) ) {
        W = W[ match(rownames(V), rownames(W)), ]
    }

    ### Functions needs:
    ## Full matrix V        mutCatalogue        components (rows) by samples (cols)    <= HAVE
    ## Left matrix W        sigCatalogue        components (rows) by signature (cols)   <= HAVE
    ## Right matrix H       expCatalogue        signature (rows) by samples (cols)     <= WANT

    # component_by_sample => NxM => N - Features, M - Samples => Component by Sample matrix
    # component_by_signature => NxL => N - Features, L - Signatures => Component by Signature matrix
    Hraw = as.matrix(LinCombDecompSigs(component_by_sample = V, component_by_signature = W))
    return(Hraw)

}

# Extracted version of core component of osCN function with fix toggle
# Loop over chromosomes to identify oscillation
testOC <- function(segTab,usefix=FALSE){

    samps = unique(segTab$sample)
    chrs = unique(segTab$chromosome)
    oscCounts = c()
    for(i in samps) {

        # Retrieve segments
        #segTab = abs_profiles[[i]]
        colnames(segTab)[4] = "segVal"

        # Loop over chromosomes to identify oscillation
        chrs = unique(segTab$chromosome)
        oscCounts = c()
        for(c in chrs) {

            currseg = as.numeric(segTab$segVal[segTab$chromosome == c])
            currseg = round(as.numeric(currseg))

            # Only take chains into consideration with a length of more than 3 elements
            if(length(currseg)>3) {
                prevval = currseg[1]
                count = 0
                for(j in 3:length(currseg)) {
                    if(currseg[j] == prevval & currseg[j] != currseg[j-1]) {
                        count = count+1
                        ## suggested fix for end of chromsome chain counts
                        ## https://github.com/markowetzlab/CINSignatureQuantification/issues/22
                        ## Included boolean to switch fix on and off
                        if(usefix){
                            if (j == length(currseg)) {
                                oscCounts = c(oscCounts, count)
                                count = 0
                            }
                        }
                    } else {
                        oscCounts = c(oscCounts,count)
                        count = 0
                    }
                    prevval = currseg[j-1]
                }
            }
        }

    }
    return(oscCounts)
}

## Report osCN chains
reportChain <- function(x,y){
    paste0(c("noFix = ","fix = "),
           c(paste0(x,collapse = ""),
             paste0(y,collapse = "")))
}

## Define function to split the samples by dCIN+ cancer type
splitSxCByGroup <- function(data=NULL,sampleCol="sample",group=NULL,groupdata=NULL,minsize=100){
    if(is.null(data)){
        stop("no data")
    }
    if(is.null(group)){
        stop("no group")
    }
    if(is.null(groupdata)){
        stop("no group data")
    }
    if(!all(c(sampleCol,group) %in% colnames(groupdata))){
        stop("missing cols to split by group")
    }

    subgroupdata <- groupdata[,c(sampleCol,group)]
    sampsplit <- split(subgroupdata,f=subgroupdata[,group])
    keepGroup <- unlist(lapply(sampsplit,FUN = function(x) nrow(x) >= minsize))
    sampsplit <- sampsplit[keepGroup]

    splitsxc <- lapply(sampsplit,FUN = function(x){
        sampleNames <- as.character(x[,sampleCol])
        subsxc <- data[rownames(data) %in% sampleNames,]
        return(subsxc)
    })


    return(splitsxc)
}

## Identify and read in files
loadNMFresults = function(PATHTOFILES, SIGMATPATTERN="W.txt", TRANSPOSED=FALSE, VERBOSE=FALSE) {

    ## Read in matrices (results from NMF runs)
    allSigFiles = list.files( path = PATHTOFILES, pattern = SIGMATPATTERN, recursive = TRUE, full.names = TRUE )
    allSamples = sapply(allSigFiles, function(x) strsplit( basename(x), "_")[[1]][1] )
    lData = lapply(allSamples, function(thisSample) {

        if(VERBOSE) { print(thisSample) }
        thisW = fread( allSigFiles[ grep(thisSample, allSigFiles) ] )
        theseRownames = thisW$V1
        thisW$V1 = NULL
        thisSigs = as.matrix(thisW)
        rownames(thisSigs) = theseRownames

        # For easier downstream analysis, transpose the matrix, so the signatures are in the rows
        if( TRANSPOSED ) {
            thisSigs = t( apply(thisSigs, 1, function(x) x/sum(x)) )
            return( thisSigs )
        } else {
            return( t(thisSigs) )
        }


    } )

    names(lData) = allSamples

    lResults = list(lData = lData, allSamples = allSamples)
    return(lResults)
}

## Determine K
detK = function(lData, OVERRIDEK=FALSE) {

    # Get number of K from all solutions
    dfKs = data.frame(sapply(lData, function(x) nrow(x)))

    ## Allow for a user-specified K
    if(isFALSE(OVERRIDEK) | ! is.numeric(OVERRIDEK)) {
        thisK = as.numeric( names(table(dfKs))[ which.max( table(dfKs) ) ] )
    } else {
        thisK = OVERRIDEK
        # Sanity check:
        if(! thisK %in% dfKs[,1]) {
            stop("Manually supplied K (via variable OVERRIDEK) is actually not present in the data.")
        }
    }

    return(thisK)
}

# From: https://stats.stackexchange.com/questions/97051/building-the-connection-between-cosine-similarity-and-correlation-in-r
cosineSimForMatrices <- function(X, corr=FALSE){
    if(corr){ X = apply(X, 2, function(x){ x-mean(x) }) }
    denom = solve(diag(sqrt(diag(t(X)%*%X))))
    return( denom%*%(t(X)%*%X)%*%denom )
}

## Identify optimal solution
idOptimalSolution = function(lData, allSamples, PATHTOFILES, LOGPATTERN="log.txt", thisK, DECISION="div", OUTPUTDIR, VERBOSE=FALSE) {

    # DECISION can be either "cost" or "div" => Result of overall cost function or Kullback-Leibler divergence

    ## Identify all log files
    allLogs = list.files( path = PATHTOFILES, pattern = LOGPATTERN, recursive = TRUE, full.names = TRUE )

    ## Loop over samples from optimal K and chose solution with best beta_divergence.
    lScores=list()
    bestCostFun = Inf
    for(thisSample in allSamples) {

        thisMat = lData[[ thisSample ]]
        if(nrow(thisMat) == thisK) {
            thisLog = fread( allLogs[ grep(pattern = thisSample, x = allLogs) ] )

            if( DECISION == "div" ) {
                thisCost = as.numeric( strsplit( thisLog[ nrow(thisLog) ]$V3, "=")[[1]][2] )
                lScores[[length(lScores)+1]] = c(thisSample, thisCost)
                #print(thisCost)
                if(thisCost < bestCostFun) {
                    if(VERBOSE) {
                        print(paste("bestSample:", thisSample,"| bestDiv:", thisCost))
                    }
                    bestSample = thisSample
                    bestCostFun = thisCost
                }
            } else {
                thisCost = as.numeric( strsplit( thisLog[ nrow(thisLog) ]$V2, "=")[[1]][2] )
                #print(thisCost)
                if(thisCost < bestCostFun) {
                    if(VERBOSE) {
                        print(paste("bestSample:", thisSample,"| bestCost:", thisCost))
                    }
                    bestSample = thisSample
                    bestCostFun = thisCost
                }
            }
        }
    }

    # For later
    dtScores=data.table(do.call(rbind, lScores))
    dtScores$V2=as.numeric(dtScores$V2)
    dtScores=dtScores[order(dtScores$V2, decreasing = FALSE),]
    colnames(dtScores) = c("Sample", as.character(DECISION))

    # Save all suitable solutions and their value
    write.table(dtScores, file.path(OUTPUTDIR, paste0("0_summary_BayesNMF_runs_all_Solutions_K", thisK, ".txt")),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

    lResults = list(dtScores = dtScores, bestSample = bestSample, thisK = thisK, decision = DECISION)
    return(lResults)
}

#### Plotting functions
## Histogram of K's
plotHistOfKs = function(lData) {
    dfKs = data.frame(sapply(lData, function(x) nrow(x)))
    colnames(dfKs) = "K"
    plotHist = ggplot(dfKs, aes(x=K)) + geom_histogram() + theme_bw() + ylab("Frequency")

    return(plotHist)
}

## Plot heatmap of cosine similarities
# Source of cosine similarity threshold of 0.85: Alexandrov et al., 2018 (bioRxiv)
plotCosineHeat = function(lData, allSamples, thisK, corThreshold = 0.85) {

    lFilt = list()
    for(thisSample in allSamples) {
        thisMat = lData[[ thisSample ]]
        if(nrow(thisMat) == thisK) {
            rownames(thisMat) = paste0( thisSample, "_", rownames(thisMat) )
            lFilt[[ thisSample ]] = thisMat
        }
    }


    matFilt = do.call(rbind, lFilt)
    corMat = cosineSimForMatrices(t(matFilt), corr = FALSE)
    corMat[ corMat < corThreshold ] = 0
    return( Heatmap(corMat, show_row_names = FALSE, show_column_names = FALSE, name = paste("K =", thisK)) )

}

makeGraphLayout = function(corMat, bestSample, thisK) {
    net = graph_from_adjacency_matrix(corMat, weighted = TRUE, mode = "undirected", diag = FALSE)

    # Get samples names
    theseSamples = sapply(rownames(corMat), function(thisSig) { strsplit(thisSig, "_")[[1]][1] })
    V(net)$Samples = theseSamples

    V(net)$color <- "#FFFFFF"
    V(net)$color[ V(net)$Samples == bestSample ] = "#228B22"

    # make TCGA larger
    V(net)$size = 3
    V(net)$size[ V(net)$Samples == bestSample ] = 10

    # label only TCGA
    V(net)$label = NA
    V(net)$label[ V(net)$Samples == bestSample ] = seq(1:thisK)

    # layout for spacing out vertices
    e = get.edgelist(net, names=FALSE)
    l = qgraph.layout.fruchtermanreingold(e, vcount = vcount(net),
                                          area=8*(vcount(net)^2),
                                          repulse.rad=(vcount(net)^3.1))

    lOut = list(net = net, l = l)
    return(lOut)
}

makeGraph = function(net, l) {
    plot(net, edge.width = 1, vertex.label.font = 2, vertex.label.color="black", layout=l)
}

## Plot hairball plot
plotHairballAndHeatmap = function(lData, allSamples, bestSample, thisK, corThreshold = 0.85) {

    # Hairball plot with highlighted best solution
    lFilt = list()
    for(thisSample in allSamples) {
        thisMat = lData[[ thisSample ]]
        if(nrow(thisMat) == thisK) {
            rownames(thisMat) = paste0( thisSample, "_", rownames(thisMat) )
            lFilt[[ thisSample ]] = thisMat
        }
    }

    matFilt = do.call(rbind, lFilt)
    corMat = cosineSimForMatrices(t(matFilt), corr = FALSE)
    corMat[ corMat < corThreshold ] = 0
    rownames(corMat) = rownames(matFilt)
    colnames(corMat) = rownames(matFilt)

    # Convert matrix to graph
    lGraph = makeGraphLayout(corMat, bestSample, thisK)

    # Plotting heatmap of best solution
    plotHeatBest = Heatmap(bestSolution, name = bestSample, cluster_rows = FALSE, cluster_columns = FALSE)

    lPlots = list(hairball = lGraph, heatmap = plotHeatBest)
    return(lPlots)
}

ViewHairball = function(lGraph) {
    makeGraph(lGraph$net, lGraph$l)
    return(lGraph)
}

## For manual solution chosing
testSolution = function(lData, allSamples, DECISION="div", bestSample, corThreshold = 0.85) {

    # Hairball plot with highlighted best solution
    lFilt = list()
    for(thisSample in allSamples) {
        thisMat = lData[[ thisSample ]]
        if(nrow(thisMat) == thisK) {
            rownames(thisMat) = paste0( thisSample, "_", rownames(thisMat) )
            lFilt[[ thisSample ]] = thisMat
        }
    }

    matFilt = do.call(rbind, lFilt)
    corMat = cosineSimForMatrices(t(matFilt), corr = FALSE)
    corMat[ corMat < corThreshold ] = 0
    rownames(corMat) = rownames(matFilt)
    colnames(corMat) = rownames(matFilt)

    # Convert matrix to graph
    lGraph = makeGraphLayout(corMat, bestSample, thisK)
    makeGraph(lGraph$net, lGraph$l)
}

## Combine plots in pdf and move files from best solution into user-selected directory
makeOutputFiles = function(OUTPUTDIR, bestSample, thisK, plotHist, plotHeat1, plotSigs, PATHTOFILES, TRANSPOSED=FALSE, SIGMATPATTERN="W.txt", EXPMATPATTERN="H.txt") {

    ## Save plots
    pdf(file = file.path(OUTPUTDIR, paste0("1_summary_BayesNMF_runs_", bestSample, "_K", thisK, ".pdf")), width = 12, height = 8)
    print(plotHist)
    print(plotHeat1)
    makeGraph(plotSigs$hairball$net, plotSigs$hairball$l)
    print(plotSigs$heatmap)
    dev.off()

    ## Identify exposure/activity and definition files for moving
    allSigFiles = list.files( path = PATHTOFILES, pattern = SIGMATPATTERN, recursive = TRUE, full.names = TRUE )
    allExpFiles = list.files( path = PATHTOFILES, pattern = EXPMATPATTERN, recursive = TRUE, full.names = TRUE )

    # Prepare exposure/activity matrix for export
    bestExpRaw = fread( allExpFiles[ grep( bestSample, allExpFiles ) ])
    expRownames = bestExpRaw$V1
    bestExpRaw$V1 = NULL
    matBestExpRaw = as.matrix(bestExpRaw)
    rownames(matBestExpRaw) = expRownames
    if( TRANSPOSED ) {
        matBest = t( apply(matBestExpRaw, 1, function(x) x/sum(x)) )
    } else {
        matBest = t( apply(matBestExpRaw, 2, function(x) x/sum(x)) )
    }


    ## Move best solution to output directory
    file.copy( allSigFiles[ grep(bestSample, allSigFiles)], file.path(OUTPUTDIR, paste0("2_Signatures_", bestSample, ".txt")))
    file.copy( allExpFiles[ grep(bestSample, allExpFiles)], file.path(OUTPUTDIR, paste0("2_Exposures_", bestSample, ".txt")))

    # Save signature matrix (normalised by default in BayesNMF algorithm)
    saveRDS(object = bestSolution, file = file.path(OUTPUTDIR, paste0("2_Signatures_", bestSample, "_normalised.rds")))

    # Save exposure matrix
    saveRDS(object = matBestExpRaw, file = file.path(OUTPUTDIR, paste0("2_Exposures_", bestSample, ".rds")))
    saveRDS(object = matBest, file = file.path(OUTPUTDIR, paste0("2_Exposures_", bestSample, "_normalised.rds")))

    return(TRUE)
}

readSigs = function(PANCANCERSIGS, name = "TCGA") {

    if(grepl(".txt", PANCANCERSIGS )) {
        pcSigs = read.table(PANCANCERSIGS, header = TRUE, row.names = 1, sep = "\t")
    } else {
        pcSigs = readRDS(PANCANCERSIGS)
    }

    if( nrow(pcSigs) > ncol(pcSigs)) pcSigs = t(pcSigs)
    rownames(pcSigs) = paste0(name, "_", rownames(pcSigs))

    return(pcSigs)

}

readExp = function(PANCANCEREXP, name = "TCGA") {

    if(grepl(".txt", PANCANCEREXP )) {
        pcExp = read.table(PANCANCEREXP, header = TRUE, row.names = 1, sep = "\t")
    } else {
        pcExp = readRDS(PANCANCEREXP)
    }

    if( nrow(pcExp) > ncol(pcExp)) pcExp = t(pcExp)
    rownames(pcExp) = paste0(name, "_", rownames(pcExp))

    return(pcExp)

}

plotHairball = function(matAll, outfile, numPCSigs, corThreshold = 0.90, cancerCols, corMethod = "pearson", plot2File = FALSE) {

    library("lsa")
    library("igraph")
    library("qgraph")
    library("RColorBrewer")

    # make network object from correlation values
    if(corMethod == "cosine") {
        # corMat = cosineSimForMatrices(t(matAll), corr = FALSE)
        corMat = cosine(t(matAll))
    } else {
        corMat = cor( t(matAll), method = corMethod )
    }
    corMat[ corMat < corThreshold ] = 0
    net = graph_from_adjacency_matrix(corMat, weighted = TRUE, mode = "undirected", diag = FALSE)

    # Add cancer type
    theseCancers = sapply(rownames(matAll), function(thisSig) { strsplit(thisSig, "_")[[1]][1] })
    V(net)$Cancer = theseCancers

    # Generate colors based on media type:
    # cols = brewer.pal(nlevels( theseCancers ), "Set2")
    cols = as.character( cancerCols$V2[ match( theseCancers, cancerCols$V1 ) ] )

    # add TCGA
    cols[ is.na(cols) ] = "#FFFFFF"
    V(net)$color <- cols

    # make TCGA larger
    V(net)$size = 3
    V(net)$size[ V(net)$Cancer == "TCGA" ] = 10

    # label only TCGA
    V(net)$label = NA
    V(net)$label[ V(net)$Cancer == "TCGA" ] = seq(1:numPCSigs)

    # layout for spacing out vertices
    e = get.edgelist(net, names=FALSE)
    l = qgraph.layout.fruchtermanreingold(e, vcount = vcount(net),
                                          area=8*(vcount(net)^2),
                                          repulse.rad=(vcount(net)^3.1))

    # prepare legend
    legendNames = unique( theseCancers )
    legendColours = as.character( cancerCols$V2[ match( unique(theseCancers), cancerCols$V1 ) ] )
    legendColours[ is.na(legendColours) ] = "#FFFFFF"

    outFinal = paste0(outfile, "-", corThreshold, ".pdf")

    if(plot2File) {
        pdf(outFinal, width = 11, height = 9)
        plot(net, edge.width = 1, vertex.label.font = 2, vertex.label.color="black", layout=l )
        legend(x=1.2, y=1.2, legend = legendNames, col = legendColours, bty = "n", pch=20 , pt.cex = 3, cex = 1, text.col="black" )
        dev.off()
    } else {
        plot(net, edge.width = 1, vertex.label.font = 2, vertex.label.color="black", layout=l )
        legend(x=1.2, y=1.2, legend = legendNames, col = legendColours, bty = "n", pch=20 , pt.cex = 3, cex = 1, text.col="black" )
    }

    ### Backup for when scaling is needed
    # l <- layout_nicely(net)
    # l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
    # plot(net, rescale=F, layout=l*1)

}

samplesPerSignature = function(lCExp,pcExp, EXPTRESH = 0.05) {

    ALLCANCERS=names(lCExp)

    pcExpNumbers = data.frame(colSums(pcExp > EXPTRESH))
    colnames(pcExpNumbers) = c("Samples")
    lExpNumbers = lapply(ALLCANCERS, function(thisCancer) {

        return( data.frame(colSums(lCExp[[thisCancer]] > EXPTRESH)) )

    } )

    csExpNumbers = do.call(rbind, lExpNumbers)
    colnames(csExpNumbers) = c("Samples")
    dtExposedSamples = rbind(pcExpNumbers, csExpNumbers)
    dtExposedSamples$Cancer = sapply(strsplit(rownames(dtExposedSamples), "_"), function(x) x[[1]])
    dtExposedSamples$Signature = sapply(strsplit(rownames(dtExposedSamples), "_"), function(x) x[[2]])

    return(dtExposedSamples)
}

removePCSimilarSigs = function(matAll, COSINETHRESH = 0.85, PREFIXPCSIGS = "TCGA") {

    # Compare and get PC sigs
    cosSigs = cosine(t(matAll))
    pcCos = cosSigs[ grepl(PREFIXPCSIGS, rownames(cosSigs)), ]

    # Identify CS sigs with no correlation to a PC sigs
    uniquePCSigs = colnames(pcCos)[ colSums(pcCos > COSINETHRESH) == 0 ]
    filtMatAll = matAll[ rownames(matAll) %in% uniquePCSigs, ]

    return(filtMatAll)
}

removeCSSimilarSigs = function(filtMatAll, dtExposedSamples, COSINETHRESH = 0.85) {

    # Prepare output
    lResults = list()

    # Compare CS sigs
    cosSigs = cosine(t(filtMatAll))

    # ID unique signatures (== 1 because of matrix diagonal)
    uniqueSigs = rownames(cosSigs)[ colSums(cosSigs > COSINETHRESH) == 1 ]
    lResults[[ length(lResults)+1 ]] = filtMatAll[rownames(filtMatAll) %in% uniqueSigs,,drop=FALSE]
    cosSigs = cosSigs[ ! rownames(cosSigs) %in% uniqueSigs, ! colnames(cosSigs) %in% uniqueSigs ]

    # Sort signatures by number of exposed samples.
    dtPresentSigs = dtExposedSamples[ rownames(dtExposedSamples) %in% rownames(cosSigs), ]
    newOrder = rownames(dtPresentSigs)[ order(dtPresentSigs$Samples, decreasing = TRUE) ]
    cosSigs = cosSigs[newOrder,newOrder]

    # Loop over signatures
    diag(cosSigs) = 0
    while(nrow(cosSigs) > 0) {

        dim(cosSigs)

        # Identify sigs
        thisName = rownames(cosSigs)[1]
        thisSig = cosSigs[1,,drop=FALSE]
        otherSigs = colnames(thisSig)[ (thisSig > COSINETHRESH) ]

        # Store winner signature and remove all from remaining matrix
        lResults[[length(lResults)+1]] = filtMatAll[rownames(filtMatAll) == thisName,,drop=FALSE]
        SOI = c(rownames(thisSig), otherSigs)
        cosSigs = cosSigs[! rownames(cosSigs) %in% SOI, ! colnames(cosSigs) %in% SOI,drop=FALSE]
    }

    csSigs = do.call(rbind, lResults)
    return(csSigs)
}

## plotCosineHeatmap function
plotCosineHeatmap <- function(x1=NULL,x2=NULL,title="",threshold=0.74,compGroup=paste0("CX",seq_len(17)),legTitle="cosine sim"){
    if(is.null(x1)){
        stop("no data")
    }
    if(is.null(x2)){
        cosinex1x2 <- cosine(t(x1))
    } else {
        cosinex1x2 <- cosine(t(rbind(x1,x2)))
    }

    cosinex1x2 <- cosinex1x2[,!colnames(cosinex1x2) %in% compGroup]
    if(!is.null(x2)){
        cosinex1x2 <- cosinex1x2[rownames(cosinex1x2) %in% compGroup,]
    }
    cosinex1x2[cosinex1x2 < threshold] = 0
    Heatmap(cosinex1x2, cluster_columns = F, cluster_rows = FALSE,
            column_title = title,name = legTitle,
            cell_fun = function(j, i, x, y, w, h, col) {
                if(round(cosinex1x2,digits = 2)[i, j] > threshold){# add text to each grid
                    grid.text(round(cosinex1x2,digits = 2)[i, j], x, y)
                }
            })
}

## getSignatureMapping function
getSignatureMapping <- function(x1=NULL,x2=NULL,threshold=0.74,compGroup=paste0("CX",seq_len(17))){
    if(is.null(x1)){
        stop("no data")
    }
    if(is.null(x2)){
        cosinex1x2 <- cosine(t(x1))
    } else {
        cosinex1x2 <- cosine(t(rbind(x1,x2)))
    }

    cosinex1x2 <- cosinex1x2[,!colnames(cosinex1x2) %in% compGroup]
    if(!is.null(x2)){
        cosinex1x2 <- cosinex1x2[rownames(cosinex1x2) %in% compGroup,]
    }
    cosinex1x2[cosinex1x2 < threshold] = 0

    maps <- do.call(rbind,apply(cosinex1x2,MARGIN = 1,FUN = function(x){
        if(any(x != 0)){
            ym <- x[which(max(x) == x)]
            nm <- names(ym)
            vm <- as.vector(ym)
            y <- x[which(x >= threshold)]
            n <- names(y)[!names(y) %in% nm]
            v <- as.vector(y[!y %in% ym])
            o <- paste0(paste0(n," = ",v,collapse = ","))
        } else {
            nm <- ""
            vm <- ""
            o <- ""
        }
        d <- data.frame(name=nm,cosineSim=vm,otherMaps=o)
        d$otherMaps[d$otherMaps == " = "] <- ""
        return(d)
    }))
    dupnames <- maps$name[which(duplicated(maps$name))]
    dupnames <- dupnames[dupnames != ""]
    maps$multimap <- ifelse(maps$name %in% dupnames,TRUE,FALSE)

    return(maps)
}

# estimateActivityThresholds
estimateThresholds <- function(feats=NULL,sigDefs=NULL,sigActs=NULL,mixtures=NULL,
                               iters=1000,SDPROP = 20,FINALWIDTH = 0.1,RANGECNAS = 0.1,
                               UNINFPRIOR = TRUE,method="drewsV2",minmaxNorm=FALSE,orderOri=FALSE,parallel=FALSE){
    if(is.null(feats)){
        stop("no data")
    }
    if(iters < 2){
        stop("require at least 2 iterations")
    }
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
            lFeatures <- CINSignatureQuantification:::simulateFeatureNoise(feats=feats,SDPROP=SDPROP,FINALWIDTH=FINALWIDTH,RANGECNAS=RANGECNAS)
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
                       acts <- calculateActivityDrewsV2(object = lSxC,SIGS = sigDefs)
                       lActs = t(apply(acts, 2, function(x) x/sum(x)))
                   })
            # progressr iter call
            p(sprintf("i=%g", i))
            list(lActs)
        }
        #doParallel::stopImplicitCluster()
    } else {
        lActsList <- list()
        for(i in 1:iters){
            message(paste0("Noise sim iteration ",i," of ",iters))
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
                       acts <- calculateActivityDrewsV2(object = lSxC,SIGS = sigDefs)
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

            dfFeat$value2[dfFeat$redraw] = rnorm(n = length(dfFeat$value[dfFeat$redraw]),
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
    dtSigs <- data.table(formatSignatureList(lSignatures),stringsAsFactors = TRUE)
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
            dtOrder = aggregate(value ~ Var1, dtCS2, median)
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
        thresh3 = qnorm(0.95, mean = modelMean, sd = sqrt(modelVar), lower.tail = TRUE)
        thresh4 = qnorm(0.99, mean = modelMean, sd = sqrt(modelVar), lower.tail = TRUE)

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
    quantile(x, probs = 0.25)
}
q75 <- function(x){
    quantile(x, probs = 0.75)
}
q95 <- function(x){
    quantile(x, probs = 0.95)
}

fitMixturePois <- function(component=NULL,seed=NULL,model_selection="BIC",min_prior=0.001,niter=1000,nrep=1,min_comp=2,max_comp=10,iters=1,cores=1){
    if(!requireNamespace("flexmix", quietly = TRUE)){
        stop("package 'flexmix' is not installed and is required to fit mixture models")
    }
    if(is.null(component)){
        stop("component is null")
    }
    if(!is.null(seed)){
        set.seed(seed)
    }
    if(min_comp >= max_comp){
        stop("minimum number of components should be less than maximum")
    }
    # capture args
    argg <- unlist(c(as.list(environment()))[-1])
    # set the column names to fix bpchrarm naming
    colnames(component) <- c("ID","value")

    control <- generateControl(min_prior = min_prior, niter = niter)

    model <- fitPoisModel(component = component$value,min_comp = min_comp,
                          max_comp = max_comp,nrep = nrep, control = control,
                          model_selection = model_selection,iters=iters,cores=cores)
    modelArg <- append(model,argg)
    return(modelArg)
}

formatPoisModel <- function(model){
    componentMeans <- modeltools::parameters(model)
    ## get prior probability for component assignment
    ## equivalent to colSums(fitted@posterior$scaled) / sum(fitted@posterior$scaled)
    componentPriors <- modeltools::prior(model)
    CompTab <- data.table::data.table(Mean = componentMeans,Weight = componentPriors)
    CompTab <- CompTab[order(CompTab$Mean)]
    return(CompTab)
}
