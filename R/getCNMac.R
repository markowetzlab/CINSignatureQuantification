getCNMac<-function(abs_profiles){
    out<-c()
    samps<-names(abs_profiles)
    for(i in samps)
    {
        # Data is pre-processed prior so should never be a QDNAseqCopyNumbers
        # when calling this function.
        # if(class(abs_profiles)=="QDNAseqCopyNumbers")
        # {
        #     segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        # }
        # else
        # {
        segTab<-abs_profiles[[i]]
        colnames(segTab)[4]<-"segVal"
        #}
        segTab$segVal[as.numeric(segTab$segVal)<0]<-0
        cn<-as.numeric(segTab$segVal)
        out<-rbind(out,cbind(ID=rep(i,length(cn)),value=cn))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}
