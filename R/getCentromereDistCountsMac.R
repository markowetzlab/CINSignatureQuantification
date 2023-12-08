getCentromereDistCountsMac<-function(abs_profiles,centromeres,chrlen){
    out<-c()
    samps<-names(abs_profiles)
    for(i in samps)
    {
        # Data is pre-processed prior so should never be a QDNAseqCopyNumbers
        # when calling this function.
        # if(class(abs_profiles)=="QDNAseqCopyNumbers")
        # {
        #     segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        # }else
        # {
        segTab<-abs_profiles[[i]]
        colnames(segTab)[4]<-"segVal"
        #}
        chrs<-unique(segTab$chromosome)
        all_dists<-c()
        for(c in chrs)
        {
            if(nrow(segTab)>1)
            {
                starts<-as.numeric(segTab$start[segTab$chromosome==c])[-1]
                segstart<-as.numeric(segTab$start[segTab$chromosome==c])[1]
                ends<-as.numeric(segTab$end[segTab$chromosome==c])
                segend<-ends[length(ends)]
                ends<-ends[-length(ends)]
                centstart<-as.numeric(centromeres[substr(centromeres[,2],4,5)==c,3])
                centend<-as.numeric(centromeres[substr(centromeres[,2],4,5)==c,4])
                chrend<-chrlen[substr(chrlen[,1],4,5)==c,2]
                ndist<-cbind(rep(NA,length(starts)),rep(NA,length(starts)))
                ndist[starts<=centstart,1]<-(centstart-starts[starts<=centstart])/(centstart-segstart)*-1
                ndist[starts>=centend,1]<-(starts[starts>=centend]-centend)/(segend-centend)
                ndist[ends<=centstart,2]<-(centstart-ends[ends<=centstart])/(centstart-segstart)*-1
                ndist[ends>=centend,2]<-(ends[ends>=centend]-centend)/(segend-centend)
		ndist<-stats::na.omit(ndist)
		ndist<-apply(ndist,1,min)

                all_dists<-rbind(all_dists,sum(ndist>0))
                all_dists<-rbind(all_dists,sum(ndist<=0))
            }
        }
        if(nrow(all_dists)>0)
        {
            out<-rbind(out,cbind(ID=i,ct1=all_dists[,1]))
        }
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}
