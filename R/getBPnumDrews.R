getBPnumDrews = function(abs_profiles,chrlen, SIZE = 10000000) {

    # Prepare looping
    out = c()
    samps = names(abs_profiles)

    # Loop over samples
    for(i in samps) {
        # Retrieve segments
        segTab = abs_profiles[[i]]
        colnames(segTab)[4] = "segVal"

        # Loop over chromosomes and identify breaks
        chrs = unique(segTab$chromosome)
        allBPnum = c()
        for(c in chrs) {
            currseg = segTab[segTab$chromosome == c,]
            intervals = seq(1, chrlen[chrlen[,1] == paste0("chr",c),2]+SIZE, SIZE)
            res = graphics::hist(as.numeric(currseg$end[-nrow(currseg)]),breaks=intervals,plot=FALSE)$counts
            allBPnum = c(allBPnum,res)
        }
        # Make sure it's really numeric
        out = rbind(out, cbind(ID = rep(i,length(allBPnum)),
                               value = as.numeric(allBPnum)))
    }

    # Prepare return
    rownames(out) = NULL
    return(data.frame(out,stringsAsFactors = FALSE))
}
