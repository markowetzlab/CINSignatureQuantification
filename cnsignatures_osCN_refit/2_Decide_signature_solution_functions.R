#### Functions

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
# or just use:
# library(lsa)
# cosine(matrix)
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


