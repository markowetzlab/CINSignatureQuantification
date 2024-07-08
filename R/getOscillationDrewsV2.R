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
