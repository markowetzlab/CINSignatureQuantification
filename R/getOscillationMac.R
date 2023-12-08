getOscillationMac<-function(abs_profiles,chrlen){
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
        oscCounts<-c()
        for(c in chrs)
        {
            currseg<-segTab$segVal[segTab$chromosome==c]
            currseg<-round(as.numeric(currseg))
            if(length(currseg)>3)
            {
                prevval<-currseg[1]
                count=0
                for(j in 3:length(currseg))
                {
                    if(currseg[j]==prevval&currseg[j]!=currseg[j-1])
                    {
                        count<-count+1
                    }else{
                        oscCounts<-c(oscCounts,count)
                        count=0
                    }
                    prevval<-currseg[j-1]
                }
            }
        }
        out<-rbind(out,cbind(ID=rep(i,length(oscCounts)),value=oscCounts))
        if(length(oscCounts)==0)
        {
            out<-rbind(out,cbind(ID=i,value=0))
        }
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}
