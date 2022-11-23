#' plotSegments
#'
#' Plot the segment data for a given sample stored in a CNQuant or SigQuant
#' class object.
#'
#' @param object CNQuant or SigQuant class object
#' @param sample vector of length 1 containing either a sample name or sample
#'   index
#' @param cn.max maximum copy number to plot - Values greater than this are
#'   truncated to the specified value
#'
#' @return plot
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
#'   segments <- plotSegments(cnobj,sample=1,cn.max=15)
#' @seealso [getSegments()]
#' @export plotSegments
#'
plotSegments <- function(object=NULL,sample=NULL,cn.max=15){
    if(is.null(object)){
        stop("No object provided, object should be a object of class CNQuant or SigQuant")
    }
    if(!class(object) %in% c("CNQuant","SigQuant")){
        stop("Object is not of class CNQuant or SigQuant")
    }
    if(is.null(sample)){
        stop("No sample specified; sample should be an integer index or name of sample contained within the provided object")
    }
    if(!is.numeric(sample) & !is.character(sample)){
        stop("Unknow sample value provided; sample should be an integer index or name of sample contained within the provided object")
    }
    samp <- getSamples(object = object)
    samp.len <- length(samp)
    if(is.numeric(sample)){
        if(sample > samp.len){
            stop(paste0("Sample index is out of bounds; Object contains ",samp.len," samples"))
        }
    }
    if(is.character(sample)){
        if(!sample %in% samp){
            stop("Sample was not found in the object provided")
        }
    }
    samp.name <- ifelse(is.numeric(sample),samp[sample],sample)
    object <- object[samp.name]
    ob.pl <- getSamplefeatures(object = object)$ploidy
    segTab <- getSegments(object = object)
    segTab$chromosome <- factor(segTab$chromosome,
                                levels = stringr::str_sort(unique(segTab$chromosome),
                                                           numeric = T))
    if(max(segTab$segVal) > cn.max){
        segTab$segVal[segTab$segVal > cn.max] <- cn.max
        ylim <- c(0,cn.max)
    } else {
        ylim <- c(0,round(max(segTab$segVal))+1)
    }
    seg.n <- nrow(segTab)
    chrom.len <- data.frame(Group.1=unique(segTab$chromosome))
    chrom.len$x.max <- stats::aggregate(segTab$end,
                                        by = list(segTab$chromosome),FUN = max)$x
    chrom.len$x.min <- stats::aggregate(segTab$start,
                                        by = list(segTab$chromosome),FUN = min)$x
    chrom.len$Group.1 <- factor(chrom.len$Group.1,
                                levels = stringr::str_sort(unique(chrom.len$Group.1),
                                                           numeric = T))
    chrom.len <- chrom.len[order(chrom.len$Group.1),]

    segTab$startf <- getcoordinates(chr = segTab$chromosome,
                                    pos = segTab$start,
                                    chrom.len = chrom.len)
    segTab$endf <- getcoordinates(chr = segTab$chromosome,
                                  pos = segTab$end,
                                  chrom.len = chrom.len)
    chrom.len$flatm <- getcoordinates(chr = chrom.len$Group.1,
                                      pos = chrom.len$x.min,
                                      chrom.len = chrom.len)
    chrom.len$flats <- getcoordinates(chr = chrom.len$Group.1,
                                      pos = chrom.len$x.max/2,
                                      chrom.len = chrom.len)
    chrom.len$flate <- getcoordinates(chr = chrom.len$Group.1,
                                      pos = chrom.len$x.max,
                                      chrom.len = chrom.len)

    title <- samp.name
    sub.title <- paste0("ploidy: ",ob.pl," | segments: ",seg.n)
    rect.col <- ifelse(seq_along(chrom.len$Group.1) %% 2 == 0,"white","grey95")

    graphics::par(mar=c(5, 4, 4, 4) + 0.2,xpd=FALSE)
    graphics::plot(NA,
         xlab="chromosome",
         ylab="absolute copy number",
         las=1,
         xlim=c(min(chrom.len$flatm),
                max(chrom.len$flate)),
         ylim=ylim,
         xaxs="i",
         xaxt="n",
         yaxp=c(ylim[1], ylim[2], ylim[2]-ylim[1]),
         yaxs="i")
    graphics::rect(xleft = chrom.len$flatm,
         xright = chrom.len$flate,
         ybottom = 0,ytop = cn.max,
         col=rect.col,
         border = NA)
    graphics::axis(1, at=chrom.len$flats, labels=chrom.len$Group.1)
    graphics::box()
    graphics::mtext(side=3, line=2, at=-0.07, adj=0, cex=1.2, title)
    graphics::mtext(side=3, line=1, at=-0.07, adj=0, cex=1, sub.title)
    graphics::abline(h = seq.int(1,cn.max-1,1),lty="dashed",col="gray50")
    graphics::segments(x0 = segTab$startf,y0 = segTab$segVal,
                       x1 = segTab$endf,y1 = segTab$segVal,lwd=3,col="blue")
}

## helper function to get coordinates from segment data
getcoordinates <- function(chr, pos, chrom.len) {
    posflat <- pos
    offset <- 0
    for (contig_ix in 1:nrow(chrom.len)) {
        on_contig <- chr == chrom.len$Group.1[contig_ix]
        posflat[on_contig] <- pos[on_contig] + offset
        offset <- offset + chrom.len$x.max[contig_ix]
    }
    posflat
}
