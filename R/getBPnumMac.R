getBPnumMac <- function(abs_profiles,chrlen){
    out <- c()
    samps <- names(abs_profiles)
    for(i in samps)
    {
        segTab <- abs_profiles[[i]]
        colnames(segTab)[4] <- "segVal"
        chrs <- unique(segTab$chromosome)
        allBPnum <- c()
        for(c in chrs)
        {
            currseg <- segTab[segTab$chromosome==c,]
            intervals <- seq(1,chrlen[chrlen[,1]==paste0("chr",c),2]+10000000,10000000)
            res <- graphics::hist(as.numeric(currseg$end[-nrow(currseg)]),breaks=intervals,plot=FALSE)$counts
            allBPnum <- c(allBPnum,res)
        }
        out <- rbind(out,cbind(ID=rep(i,length(allBPnum)),value=allBPnum))
    }
    rownames(out) <- NULL
    data.frame(out,stringsAsFactors = F)
}
