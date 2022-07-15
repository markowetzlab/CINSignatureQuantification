getCancerSpecificSignatures <- function(x){
    cancerSpecificSignatures = get(utils::data("cancerSpecificSignatures",envir = environment()),)
    if(!x %in% names(cancerSpecificSignatures)){
        stop("Unknown cancer subset '",x,"' provided. The following have cancer-specific signature groups:\n",
             paste0(names(cancerSpecificSignatures),collapse=", "),call. = F)
    } else {
        sigs <- as.character(cancerSpecificSignatures[[x]])
        return(sigs)
    }

}
