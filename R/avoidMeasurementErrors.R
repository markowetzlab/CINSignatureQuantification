avoidMeasurementErrors = function(dtSmooth) {
    dtSmooth$segVal[ dtSmooth$segVal < 0 ] = 0
    dtSmooth$sample = factor(dtSmooth$sample)
    return(dtSmooth)
}
