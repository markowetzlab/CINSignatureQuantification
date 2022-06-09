removeQuietSamples = function(dtSmooth, DCIN = 20) {

    # Identify CNAs per sample
    dtCNAs = dtSmooth[ dtSmooth$segVal != 2, ]
    quietSamples = names(table(dtCNAs$sample))[ table(dtCNAs$sample) < DCIN ]

    # Remove quiet samples
    dtSmooth = dtSmooth[ ! dtSmooth$sample %in% quietSamples, ]
    dtSmooth$sample = factor(dtSmooth$sample)

    return(dtSmooth)
}
