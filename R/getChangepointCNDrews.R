getChangepointCNDrews = function(abs_profiles, allowedError = 0.1, rmNorm = FALSE) {

    # Prepare looping
    out = c()
    samps = names(abs_profiles)
    # Loop over samples
    for(i in samps) {

        # Retrieve segments
        segTab = abs_profiles[[i]]
        colnames(segTab)[4] = "segVal"

        # Initiate and prepare looping over chromosomes
        segTab$segVal = as.numeric(segTab$segVal)
        segTab$segVal[segTab$segVal<0] = 0
        chrs = unique(segTab$chromosome)
        allcp = c()

        # Loop over chromosomes
        for(c in chrs) {
            currseg = as.numeric(segTab$segVal[segTab$chromosome==c])
            firstSeg = abs(2 - currseg[1] )
            # As we look only at the left end of a CNA, we might miss a changepoint at the beginning of the p-arm
            # That's why we check manually but only regard this value if it is higher than an allowed error rate.
            if(firstSeg <= allowedError) {
                theseChanges = abs(currseg[-1]-currseg[-length(currseg)])
                if(rmNorm) { theseChanges = theseChanges[ currseg[-1] != 2 ] }
                allcp = c(allcp, theseChanges)
            } else {
                theseChanges = c( firstSeg, abs(currseg[-1]-currseg[-length(currseg)]) )
                if(rmNorm) { theseChanges = theseChanges[ currseg != 2 ] }
                allcp = c(allcp, theseChanges)
            }

        }
        if(length(allcp)==0) {
            allcp = 0 #if there are no changepoints
        }
        # Make sure it's really numeric
        out = rbind(out,cbind(ID=rep(i,length(allcp)),value=as.numeric(allcp)))
    }

    # Prepare return
    rownames(out) = NULL
    return(data.frame(out,stringsAsFactors = FALSE))
}
