getSegsizeMac<-function(abs_profiles){
    out <- c()
    samps <- names(abs_profiles)
    for(i in samps)
    {
        # Data is pre-processed prior so should never be a QDNAseqCopyNumbers
        # when calling this function.
        # if(class(abs_profiles)=="QDNAseqCopyNumbers")
        # {
        #     segTab <- getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        # }
        # else
        # {
        segTab <- abs_profiles[[i]]
        colnames(segTab)[4] <- "segVal"
        #}
        segTab$segVal[as.numeric(segTab$segVal)<0]<-0
        seglen <- (as.numeric(segTab$end)-as.numeric(segTab$start))
        seglen <- seglen[seglen>0]
        out <- rbind(out,cbind(ID=rep(i,length(seglen)),value=seglen))
    }
    rownames(out) <- NULL
    data.frame(out,stringsAsFactors = F)
}
