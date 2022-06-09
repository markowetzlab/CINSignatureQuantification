scaleByModel = function(H, lModel) {
    # First go over columns to subtract mean, then go over columns to divide by scale
    Hscaled = sweep(
                sweep(H, 2, lModel$mean, FUN = '-'),
                2, lModel$scale, FUN = "/")

    return(Hscaled)
}
