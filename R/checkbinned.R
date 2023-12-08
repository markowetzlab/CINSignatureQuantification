checkbinned <- function(segTable){
    t.chr <- segTable$chromosome[1]
    t.end <- segTable$end[segTable$chromosome == t.chr]
    t.start <- segTable$start[segTable$chromosome == t.chr]
    startend.len <- length(unique(t.end - t.start))
    if(startend.len < 2){
        return(TRUE)
    } else {
        return(FALSE)
    }
}
