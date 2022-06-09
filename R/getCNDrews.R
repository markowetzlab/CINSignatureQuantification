getCNDrews = function(abs_profiles, rmNorm = FALSE) {

    # Prepare looping
    out = c()
    samps = names(abs_profiles)
    # Loop over samples
    for(i in samps) {

        # Retrieve segments
        segTab = abs_profiles[[i]]
        colnames(segTab)[4] = "segVal"

        segTab$segVal[as.numeric(segTab$segVal)<0] = 0
        # If wished, don't consider normal segments
        if(rmNorm) { segTab = segTab[ segTab$segVal != 2, ] }
        cn = as.numeric(segTab$segVal)
        # Make sure it's really numeric.
        out = rbind(out,cbind(ID=rep(i,length(cn)),value=as.numeric(cn)))
    }
    # Prepare return
    rownames(out) = NULL
    return(data.frame(out,stringsAsFactors = FALSE))
}
