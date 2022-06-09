getSegsizeDrews = function(abs_profiles, rmNorm = FALSE) {

    # Prepare looping
    out = c()
    samps = names(abs_profiles)

    # Loop over samples
    for(i in samps) {

        # Get segments
        segTab = abs_profiles[[i]]
        colnames(segTab)[4] = "segVal"

        # Make sure segment values are numeric
        segTab$segVal = as.numeric(segTab$segVal)

        # If wished, don't consider normal segments
        if(rmNorm) { segTab = segTab[ segTab$segVal != 2, ] }

        # Avoiding potential artefact
        segTab$segVal[segTab$segVal<0] = 0
        seglen = segTab$end-segTab$start
        seglen = seglen[seglen>0]

        # Double tap.
        out = rbind(out,cbind(ID=rep(i,length(seglen)),value=as.numeric(seglen)))
    }

    # Prepare return
    rownames(out) = NULL
    return(data.frame(out,stringsAsFactors = FALSE))

}
