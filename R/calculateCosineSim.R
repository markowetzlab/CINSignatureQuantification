#' calculateCosineSim
#'
#' MAKE into METHOD
#'
#' Function to compute cosine similarity between various input types and
#' additionally perform permutation testing for cosine similarity significance
#' testing.
#'
#' @param x A matrix or vector (must be matrix if y is NULL)
#' @param y Optional matrix or vector (default: NULL)
#' @param perm Logical to enable permutation testing (default: FALSE)
#' @param nperms Number of permutations to compute (default: 1000)
#' @param alternative tailed test type (default: "two.sided"). Additionally can
#'   run either "greater" or "less" for corresponding one-sided tests
#'
#' @returns Either a vector, matrix, or list depending on input and 'perm' option
#' @export
#'
#' @examples
calculateCosineSim <- function(x=NULL,y=NULL,perm=FALSE,nperms=1000,alternative="two.sided"){

    if(is.null(x)){
        stop("no x specified")
    }

    if(!is.numeric(x) & !is.matrix(x)){
        stop("x must be numeric or matrix")
    }

    if(is.vector(x) & is.matrix(y)){
        stop("x must be matrix if y is null")
    }

    if(!is.null(y)){
        if(!is.numeric(y) & !is.matrix(y)){
            stop("y must be numeric or matrix")
        }
    }

    ## arg checking
    rlang::arg_match(alternative,c("two.sided","greater","less"),multiple = FALSE)
    stopifnot(is.logical(perm))
    stopifnot(is.numeric(nperms))
    stopifnot(all(nperms == floor(nperms)))

    # Set method based on input type
    assign("method",setCosineSimMethod(x,y))
    #print(method)
    switch(method,
           mat={
               cosim <- lsa::cosine(x)
               if(perm){
                   y <- x
                   co <- runMatrixPermutations(x,y,nperms=nperms,alternative=alternative)
               }
           },
           matmat={
               cosim <- apply(x,MARGIN = 2,FUN = function(z) lsa::cosine(z,y))
               if(perm){
                   co <- runMatrixPermutations(x,y,nperms=nperms,alternative=alternative)
               }
           },
           matvec={
               cosim <- lsa::cosine(x = y,y = x)
               if(perm){
                   co <- apply(x,MARGIN = 2,
                               FUN = function(z) cosinePermTest(z,y,
                                                                nperms=nperms,
                                                                alternative=alternative)["p.value"])
               }
           },
           vecmat={
               c <- colnames(y)
               cosim <- lsa::cosine(x,y)
               names(cosim) <- c
               if(perm){
                   co <- apply(y,MARGIN = 2,
                               FUN = function(z) cosinePermTest(z,x,
                                                                nperms=nperms,
                                                                alternative=alternative)["p.value"])
               }
           },
           vecvec={
               cosim <- as.numeric(lsa::cosine(x,y))
               co <- cosinePermTest(x,y,
                                    nperms=nperms,
                                    alternative=alternative)["p.value"]
           })
    if(perm){
        dimnames(co) <- dimnames(cosim)
        return(list(cosine=cosim,
                    signif=co,
                    params=c(alternative=alternative,
                             nperms=nperms,method=method)))
    } else {
        return(cosim)
    }
}

cosinePermTest <- function(x,y,nperms=1000,alternative="two.sided"){
    # based on coreGx cosinePerm function - GPLv3 license
    # https://github.com/bhklab/CoreGx/blob/main/R/cosinePerm.R
    rlang::arg_match(alternative,c("two.sided","greater","less"),multiple = FALSE)

    stopifnot(is.numeric(nperms) )
    stopifnot(all(nperms == floor(nperms)))

    if(nperms < 10){
        stop("number of permutations should be greater than 10")
    }
    res <- c(cosine=lsa::cosine(x,y),p.value=NA)

    ss <- sapply(seq_len(nperms),FUN = function(x,a,b){
        a <- sample(a)
        b <- sample(b)
        s <- lsa::cosine(a,b)
        return(s)
    }, a = x, b = y)

    lesser <- sum(ss < res["cosine"])
    greater <- sum(ss > res["cosine"])
    len <- sum(!is.na(ss))

    switch(alternative, two.sided = {
        res["p.value"] <- 2 * (min(lesser,greater) / len)
    }, less = {
        res["p.value"] <- lesser / len
    }, greater = {
        res["p.value"] <- greater / len
    })
    # Adjust for minimal p-value due to permutation limit
    if(!is.na(res["p.value"])){
        if (res["p.value"] == 0) {
            res["p.value"] <- 1 / (nperms + 1)
        }
    }
    return(res)
}

setCosineSimMethod <- function(x,y){
    if(is.matrix(x) && is.null(y)){
        method <- "mat"
    } else if(is.matrix(x) && is.matrix(y)){
        method <- "matmat"
    } else if(is.vector(y) && is.matrix(x)){
        method <- "matvec"
    } else if(is.vector(x) && is.matrix(y)){
        method <- "vecmat"
    } else if(is.vector(x) && is.vector(y)){
        method <- "vecvec"
    }
    return(method)
}

runMatrixPermutations <- function(x,y,nperms=1000,alternative="two.sided"){
    co <- array(0, c(ncol(y), ncol(x)))
    for(i in 1:ncol(y)){
        for(j in 1:ncol(x)){
            co[i, j] = cosinePermTest(y[, i], x[, j],
                                      nperms=nperms,
                                      alternative=alternative)["p.value"]
        }
    }
    return(co)
}
