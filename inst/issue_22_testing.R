# test os chain bug
## https://github.com/markowetzlab/CINSignatureQuantification/issues/22

segTab_internal_chains <- data.frame(chromosome=c(rep("chr1",times=7),rep("chr2",times=7)),
                 start=rep(c(1,100,200,300,400,500,600),times=2),
                 end=rep(c(99,199,299,399,499,599,699),times=2),
                 segVal=c(c(2,2,4,3,4,2,2),c(2,2,4,3,4,2,2)),
                 sample=c(rep("sample1",times=14)))

segTab_chr2_end_chains <- data.frame(chromosome=c(rep("chr1",times=7),rep("chr2",times=7)),
                                     start=rep(c(1,100,200,300,400,500,600),times=2),
                                     end=rep(c(99,199,299,399,499,599,699),times=2),
                                     segVal=c(c(2,2,4,3,4,2,2),c(2,2,4,3,4,3,4)),
                                     sample=c(rep("sample1",times=14)))

segTab_twochains_chains <- data.frame(chromosome=c(rep("chr1",times=7),rep("chr2",times=7)),
                                     start=rep(c(1,100,200,300,400,500,600),times=2),
                                     end=rep(c(99,199,299,399,499,599,699),times=2),
                                     segVal=c(c(4,3,4,2,3,4,3),c(4,3,4,3,4,3,4)),
                                     sample=c(rep("sample1",times=14)))

# Loop over chromosomes to identify oscillation
testOC <- function(segTab,usefix=FALSE){

    samps = unique(segTab$sample)
    chrs = unique(segTab$chromosome)
    oscCounts = c()
    for(i in samps) {

        # Retrieve segments
        #segTab = abs_profiles[[i]]
        colnames(segTab)[4] = "segVal"

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
                        ## suggested fix for end of chromsome chain counts
                        ## https://github.com/markowetzlab/CINSignatureQuantification/issues/22
                        if(usefix){
                            if (j == length(currseg)) {
                                oscCounts = c(oscCounts, count)
                                count = 0
                            }
                        }
                    } else {
                        oscCounts = c(oscCounts,count)
                        count = 0
                    }
                    prevval = currseg[j-1]
                }
            }
        }

    }
    return(oscCounts)
}

int_chain_nofix <- testOC(segTab = segTab_internal_chains)
int_chains_fix <- testOC(segTab = segTab_internal_chains,usefix = T)
plotSegments(createCNQuant(segTab_internal_chains),sample = 1)

chr2_end_chain_nofix <- testOC(segTab = segTab_chr2_end_chains)
chr2_end_chains_fix <- testOC(segTab = segTab_chr2_end_chains,usefix = T)
plotSegments(createCNQuant(segTab_chr2_end_chains),sample = 1)

twochains_chain_nofix <- testOC(segTab = segTab_twochains_chains)
twochains_chains_fix <- testOC(segTab = segTab_twochains_chains,usefix = T)
plotSegments(createCNQuant(segTab_twochains_chains),sample = 1)

segTab_twochains_chains

library(CINSignatureQuantification)
library(ggplot2)

q_noFix <- quantifyCNSignatures(TCGA_478_Samples_SNP6_GOLD)
# manually load fixed version here
q_Fix <- quantifyCNSignatures(TCGA_478_Samples_SNP6_GOLD)


qFeatsNofix <- getFeatures(q_noFix)$os
qFeatsfix <- getFeatures(q_Fix)$os
ks.test(qFeatsNofix$value,qFeatsfix$value)

sum(qFeatsNofix$value > 0)
sum(qFeatsfix$value > 0)

chain_counts <- data.frame(getOS_method=factor(c("original","fixed"),levels = c("original","fixed")),
                           chains=c(sum(qFeatsNofix$value > 0),sum(qFeatsfix$value > 0)))

ggplot(chain_counts) +
    geom_col(aes(getOS_method,chains,fill=getOS_method)) +
    geom_text(aes(getOS_method,chains-400,
                  label=paste0(chains," (",(round(chains / sum(qFeatsNofix$value > 0),digits = 3))*100,"%)"))) +
    labs(title = "Chain counts (TCGA 478 gold standard)") +
    theme_bw() +
    theme(legend.position = "bottom",axis.title.x = element_blank())

chain_lengths <- data.frame(getOS_method=factor(c("original","fixed"),levels = c("original","fixed")),
                           chain_lengths=c(sum(qFeatsNofix$value),sum(qFeatsfix$value)))

ggplot(chain_lengths) +
    geom_col(aes(getOS_method,chain_lengths,fill=getOS_method)) +
    geom_text(aes(getOS_method,chain_lengths-400,
                  label=paste0(chain_lengths," (",(round(chain_lengths / sum(qFeatsNofix$value),digits = 3))*100,"%)"))) +
    labs(title = "Chain lengths (TCGA 478 gold standard)") +
    theme_bw() +
    theme(legend.position = "bottom",axis.title.x = element_blank())


chain_density <- data.frame(getOS_method=factor(c(rep("original",times=length(qFeatsNofix$value)),
                                                 rep("fixed",times=length(qFeatsfix$value))),c("original","fixed")),
                                               value=c(qFeatsNofix$value,qFeatsfix$value))

#chain_density <- chain_density[chain_density$value != 0,]

summary(chain_density$value[chain_density$getOS_method == "original"])
summary(chain_density$value[chain_density$getOS_method == "fixed"])

ggplot(chain_density[chain_density$value < 10,]) +
    geom_density(aes(value,fill=getOS_method),alpha=0.4) +
    labs(title = "Chain length distribution (TCGA 478 gold standard)") +
    theme_bw() +
    theme(legend.position = "bottom",axis.title.x = element_blank()) +
    facet_wrap(.~getOS_method)


qDiff <- getActivities(q_Fix) - getActivities(q_noFix)
barplot(colMeans(qDiff))
barplot(colMeans(qDiff) / colMeans(getActivities(q_noFix)))
