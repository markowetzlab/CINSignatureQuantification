checkSegValRounding <- function(x){
    x <- as.numeric(x)
    if(all(x == round(x,digits = 0))){
        return(TRUE)
    } else {
        return(FALSE)
    }
}
