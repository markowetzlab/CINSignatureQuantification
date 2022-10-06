getCentromereDistCountsDrews = function(abs_profiles,centromeres,chrlen) {

    # Prepare looping
    out = c()
    samps = names(abs_profiles)
    # Loop over samples
    for(i in samps) {

        # Retrieve segments
        segTab = abs_profiles[[i]]
        colnames(segTab)[4] = "segVal"

        # Loop over chromosomes
        chrs = unique(segTab$chromosome)
        all_dists = c()
        for(c in chrs) {
            if(nrow(segTab) > 1) {
                starts = as.numeric(segTab$start[segTab$chromosome==c])[-1]
                segstart = as.numeric(segTab$start[segTab$chromosome==c])[1]
                ends = as.numeric(segTab$end[segTab$chromosome==c])
                segend = ends[length(ends)]
                ends = ends[-length(ends)]
                centstart = as.numeric(centromeres[substr(centromeres[,2],4,5)==c,3])
                centend = as.numeric(centromeres[substr(centromeres[,2],4,5)==c,4])
                chrend = chrlen[substr(chrlen[,1],4,5)==c,2]
                ndist = cbind(rep(NA,length(starts)),rep(NA,length(starts)))
                ndist[starts<=centstart,1] = (centstart-starts[starts<=centstart])/(centstart-segstart)*-1
                ndist[starts>=centend,1] = (starts[starts>=centend]-centend)/(segend-centend)
                ndist[ends<=centstart,2] = (centstart-ends[ends<=centstart])/(centstart-segstart)*-1
                ndist[ends>=centend,2] = (ends[ends>=centend]-centend)/(segend-centend)
		ndist = na.omit(ndist)
		ndist = apply(ndist,1,min)

                all_dists = rbind(all_dists,sum(ndist>0))
                all_dists = rbind(all_dists,sum(ndist<=0))
            }
        }
        if(nrow(all_dists)>0) {
            # Make sure it's really numeric
            out = rbind(out,cbind(ID=i,ct1=as.numeric(all_dists[,1])))
        }
    }

    # Prepare return
    rownames(out) = NULL
    return(data.frame(out,stringsAsFactors = FALSE))
}
