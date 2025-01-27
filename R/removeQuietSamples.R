removeQuietSamples = function(dtSmooth, DCIN = 20L) {
    if(!is.numeric(DCIN)){
        stop(paste0("Provided CN event filter non-numeric"))
    } else if(DCIN < 0) {
        stop(paste0("Provided CN event filter needs to be positive"))
    } else if(DCIN %% 1 != 0){
        stop(paste0("Provided CN event filter can not be a float"))
    }
    # Identify CNAs per sample
    dtCNAs = dtSmooth[ dtSmooth$segVal != 2, ]
    quietSamples = names(table(dtCNAs$sample))[ table(dtCNAs$sample) < DCIN ]

    # Remove quiet samples
    dtSmooth = dtSmooth[ ! dtSmooth$sample %in% quietSamples, ]
    dtSmooth$sample = factor(dtSmooth$sample)

    return(dtSmooth)
}
