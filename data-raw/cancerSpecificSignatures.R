library(lsa)
library(stringr)

get_reference_sigs <- function(denovo=NULL){
    if(is.null(denovo)){
        stop("no denovo set of signatures provided")
    }
    #get(data("Drews2022_TCGA_Signatures"))
    load("data/Drews2022_TCGA_Signatures.rda")
    Drews2022_TCGA_Signatures.t <- t(Drews2022_TCGA_Signatures)

    cosine_mat <- apply(denovo,MARGIN = 2,FUN = function(x){
        apply(Drews2022_TCGA_Signatures.t,MARGIN = 2,FUN = function(y){
            lsa::cosine(x,y)
        })
    })

    #print(cosine_mat)

    cosine_mat[cosine_mat < 0.74] <- NA

    rowname_idx <- apply(cosine_mat,MARGIN = 2,FUN = function(z) ifelse(all(is.na(z)),NA,which(max(z,na.rm = T) == z)))

    mapped_sigs <- rownames(cosine_mat)[rowname_idx]

    colnames(denovo) <- mapped_sigs

    if(any(is.na(mapped_sigs))){
        new_names <- names(which(is.na(rowname_idx)))
        colnames(denovo)[which(is.na(rowname_idx))] <- new_names
    }

    return(denovo)
}


# css.path <- list.files(path = "inst/",pattern = "4_Signature*",full.names = T)
# css.names <- gsub(x = css.path,pattern = "inst/4_Signatures_(.+)\\.txt",replacement = "\\1")
#
# css.list <- do.call(list,lapply(css.path,FUN = function(x){
#     y <- as.matrix(read.table(x,header = T,sep = "\t",row.names = 1))
#     return(y)
# }))
# names(css.list) <- css.names

css.list <- readRDS("inst/ccs.list.RDS")

mapped_sigs <- lapply(css.list,FUN = function(x){
    sig.names <- unique(colnames(get_reference_sigs(denovo = x)))
    return(str_sort(sig.names,numeric=T))
})

all(paste0("CX",1:17) %in% unlist(mapped_sigs)) ## All present
any(!unlist(mapped_sigs) %in% paste0("CX",1:17)) ## No extra signatures

cancerSpecificSignatures = mapped_sigs

rm(list=ls()[ls() != "cancerSpecificSignatures"])

usethis::use_data(cancerSpecificSignatures,overwrite = TRUE)

