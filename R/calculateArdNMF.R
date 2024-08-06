#' calculateArdNMF
#'
#' @param object A object of class CNQuant with Sample-by-component matrix data
#' @param runs Number of runs to perform from which to select best ARD-NMF
#'   solution (Default: 1)
#' @param max_iter Maximum number of iterations for each run to converge
#'   (Default: 10000)
#' @param regularisation_W Regularisation term for W provided to ARD-NMF
#'   (Default: "L1")
#' @param regularisation_H Regularisation term for W provided to ARD-NMF
#'   (Default: "L1")
#' @param additional.args A list of additional arguments to pass to ARD-NMF
#'   which are checked against all available arguments
#' @param parallel Use doFuture parallelisation to implement runs set by future::plan()
#'
#' @details This function computes ARD-NMF on the provided sample-by-component
#'   matrix using a GPU-accelerated implementation originally in the
#'   \link[=https://github.com/broadinstitute/SignatureAnalyzer-GPU]{SignatureAnalyzer-GPU}
#'   repository. The python code is implemented via the
#'   \link[=https://cloud.r-project.org/web/packages/reticulate/index.html]{reticulate}
#'   package which set handles virtual environment set up, argument passing, and
#'   returning ARD-NMF outputs.
#'
#'   parallelisation is implemented via the doFuture package, then plan() for which can be
#'   specified in the global environment prior to running the function.
#'   Additionally, function calls are included to allow for the use of the
#'   progressr package to provide progress bars for both sequential and parallel
#'   processing.
#'
#'
#' @return list of ARD-NMF solutions equal to the number of runs
#' @export calculateArdNMF
#'
calculateArdNMF <- function(data=NULL,runs=1,max_iter=10000,regularisation_W="L1",
                            regularisation_H="L1",additional.args=NULL,parallel=FALSE){
    ## reticulate implementation of ardnmf
    # conda env
    # pytorch torchvision torchaudio pytorch-cuda=11.8 numpy pandas pyarrow
    # scikit-image scikit-learn scipy feather-format -c pytorch -c nvidia
    # Code here uses modified version of SignatureAnalyzer-GPU with modification to return data to R session rather than a flat file
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop(
            "Package \"reticulate\" must be installed to use ARD-NMF provided by SignatureAnalyzer-GPU.",
            call. = FALSE
        )
    }
    if(is.null(data)){
        stop("data not provided")
    } else {
        data <- as.data.frame(t(data))
    }
    # Function to check and setup conda environment - needs more testing
    ## RETICULATE version 1.38
    install_ardnmf()
    # Output directory is not used as object is returned to R but script still
    # requires an argument to be provided so a random 6 character string is used.
    randomOut <- paste0(sample(c(LETTERS,c(0:9)),size = 6,replace = T),collapse = "")
    cmdargs <- list(max_iter=max_iter,output_dir=randomOut,prior_on_W=regularisation_W,prior_on_H=regularisation_H)

    if(!is.null(additional.args)){
        if(is.list(additional.args)){
            cmdargs <- append(cmdargs,additional.args)
        } else {
            stop("Additional options should be provided as a list - e.g 'additional.args = list(phi=1,a=0.5)'")
        }
    }

    nmf_options <- parse_ard_nmf_args(args = cmdargs)

    if(runs > 1){
        if(parallel) {
            if (!requireNamespace("doFuture", quietly = TRUE)) {
                stop(
                "Package \"doFuture\" must be installed to use multiple threads/cores.",
                call. = FALSE)
            }

            message(paste0("running doFuture foreach for ",runs," runs over ",future::nbrOfWorkers(), " threads"))
            # Multi-thread usage
            `%dofuture%` <- doFuture::`%dofuture%`
            # provide progressr call if wanted
            p <- progressr::progressor(steps = runs)
            ard_runs <- foreach::foreach(r=1:runs,.options.future = list(seed = TRUE,packages = c("CINSignatureQuantification"))) %dofuture% {
                ard_run <- run_ard_nmf(data = data,nmf_options = nmf_options)
                p(sprintf("i=%g\n", r))
                list(ard_run)
            }
        } else {
            message(paste0("running sequentially for ",runs))
            ard_runs <- list()
            p <- progressr::progressor(steps = runs)
            for(r in seq_len(runs)){
                cat(paste0("Run ",r," of ",runs,"\n"))
                ard_run <- run_ard_nmf(data = data,nmf_options = nmf_options)
                p(sprintf("i=%g\n", r))
                ard_runs <- append(ard_runs,list(ard_run))
            }
        }
    } else {
        ard_runs <- list(run_ard_nmf(data = data,nmf_options = nmf_options))
    }
    names(ard_runs) <- paste0("run",seq_len(length(ard_runs)))
    return(ard_runs)
}

parse_ard_nmf_args <- function(args){
    allOpts <- c("K0","max_iter","del_","tolerance","phi","a","b","objective","prior_on_W","prior_on_H",
                 "output_dir","labeled","report_frequency","dtype","force_use_val_set","force_no_val_set")
    if(!all(names(args) %in% allOpts)){
        stop("unknown or supported args provided to ARD-NMF")
    }
    args <- unlist(c(rbind(paste0("--",names(args)),args)))
    #print(args)
    if(length(args) %% 2 != 0){
        stop(paste0("unpaired args given in",args))
    }

    args <- c(args,"--labeled")
    args[args != ""]
    return(args)
}

run_ard_nmf <- function(data=NULL,nmf_options=NULL){
    ## source ARD NMF python functions via reticulate
    reticulate::source_python(system.file("ardnmf","SignatureAnalyzer-GPU.py", package = "CINSignatureQuantification"))
    ## Run the main() function with minor alterations to return object to R rather than write to disk
    tempout <- main(data = reticulate::r_to_py(data),arglist = reticulate::r_to_py(nmf_options))
    ## Format returned nested list into single list
    ardNMFout <- format_ard_nmf(tempout)
    return(ardNMFout)
}

format_ard_nmf <- function(x){
    ## Extract first 3 elements - nsgis,W and H mat and append iter report
    y <- x[[1]]

    y$report <- as.data.frame(apply(x[[2]],MARGIN = 2,
                                    FUN = function(x) as.numeric(unlist(x))))
    names(y) <- c("nsigs","W_df","H_df","report")
    return(y)
}

install_ardnmf <- function(envname = "r-ardnmf",
                           new_env = identical(envname, "r-ardnmf"),
                           forceInstall=FALSE){

    if(forceInstall){
        if(new_env && reticulate::condaenv_exists(envname)){
            reticulate::conda_remove(envname)
        }
    }
    if(!reticulate::condaenv_exists(envname)){
        conda.packages <- c("pytorch","torchvision","torchaudio","pytorch-cuda=11.8",
                            "numpy","pandas","pyarrow","scikit-image","scikit-learn",
                            "scipy","feather-format")
        channels <- c("conda-forge","pytorch","nvidia")
        reticulate::conda_create(envname = envname,packages = conda.packages,channel = channels,additional_create_args = "--quiet")
    }

    reticulate::use_condaenv(envname)
}
