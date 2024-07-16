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

# # Replaced with "cosine" from package "lsa"
# # From: https://stats.stackexchange.com/questions/97051/building-the-connection-between-cosine-similarity-and-correlation-in-r
# cosineSimForMatrices <- function(X, corr=FALSE){
#     if(corr){ X = apply(X, 2, function(x){ x-mean(x) }) }
#     denom = solve(diag(sqrt(diag(t(X)%*%X))))
#     return( denom%*%(t(X)%*%X)%*%denom )
# }

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


samplesPerSignature = function(lCExp, EXPTRESH = 0.05) {

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

