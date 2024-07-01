# test os chain bug

segTab <- data.frame(chromosome=c(rep("chr1",times=6),rep("chr2",times=6)),
                 start=rep(c(1,100,200,300,400,500),times=2),
                 end=rep(c(99,199,299,399,499,599),times=2),
                 segVal=c(c(2,2,4,3,4,2),c(2,2,2,3,4,3)),
                 sample=c(rep("sample1",times=12)))
segTab

# Loop over chromosomes to identify oscillation
chrs = unique(segTab$chromosome)
oscCounts = c()
for(c in chrs) {

    currseg = as.numeric(segTab$segVal[segTab$chromosome == c])
    currseg = round(as.numeric(currseg))

        # Only take chains into consideration with a length of more than 3 elements
    if(length(currseg)>3) {
        prevval = currseg[1]
        count = 0
        for(j in 3:length(currseg)) {
            if(currseg[j] == prevval & currseg[j] != currseg[j-1]) {
                count = count+1
                # if (j == length(currseg)) {
                #     oscCounts = c(oscCounts, count)
                #     count = 0
                # }
            } else {
                oscCounts = c(oscCounts,count)
                count = 0
            }
            prevval = currseg[j-1]
        }
    }
}



qDiff <- abs(getActivities(q_fix) - getActivities(q_nofix))
