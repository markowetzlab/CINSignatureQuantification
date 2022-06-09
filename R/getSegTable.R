getSegTable<-function(x){
    if(class(x)=="QDNAseqCopyNumbers"){
        sn<-Biobase::assayDataElement(x,"segmented")
        fd <- Biobase::fData(x)
        fd$use -> use
        fdfiltfull<-fd[use,]
        sn<-sn[use,]
        segTable<-c()
        for(s in colnames(sn)){
            for(c in unique(fdfiltfull$chromosome))
            {
                snfilt<-sn[fdfiltfull$chromosome==c,colnames(sn) == s]
                fdfilt<-fdfiltfull[fdfiltfull$chromosome==c,]
                sn.rle<-rle(snfilt)
                starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
                ends <- cumsum(sn.rle$lengths)
                lapply(1:length(sn.rle$lengths), function(s) {
                    from <- fdfilt$start[starts[s]]
                    to <- fdfilt$end[ends[s]]
                    segValue <- sn.rle$value[s]
                    c(fdfilt$chromosome[starts[s]], from, to, segValue)
                }) -> segtmp
                segTableRaw <- data.frame(matrix(unlist(segtmp), ncol=4, byrow=T),sample = rep(s,times=nrow(matrix(unlist(segtmp), ncol=4, byrow=T))),stringsAsFactors=F)
                segTable<-rbind(segTable,segTableRaw)
            }
        }
        colnames(segTable) <- c("chromosome", "start", "end", "segVal","sample")
        return(segTable)
    } else {
        # NON QDNASEQ BINNED DATA FUNCTION

    }
}
