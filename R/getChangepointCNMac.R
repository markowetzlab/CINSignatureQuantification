getChangepointCNMac<-function(abs_profiles){
    out<-c()
    samps<-names(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }
        else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        segTab$segVal[as.numeric(segTab$segVal)<0]<-0
        chrs<-unique(segTab$chromosome)
        allcp<-c()
        for(c in chrs)
        {
            currseg<-as.numeric(segTab$segVal[segTab$chromosome==c])
            allcp<-c(allcp,abs(currseg[-1]-currseg[-length(currseg)]))
        }
        if(length(allcp)==0)
        {
            allcp<-0 #if there are no changepoints
        }
        out<-rbind(out,cbind(ID=rep(i,length(allcp)),value=allcp))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}
