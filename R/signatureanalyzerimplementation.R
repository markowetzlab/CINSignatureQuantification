## reticulate implementation of ardnmf
# used conda env
# micromamba create name bayesnmf2
    # pytorch torchvision torchaudio pytorch-cuda=11.8 numpy pandas pyarrow
    # scikit-image scikit-learn scipy feather-format
    # -c pytorch -c nvidia
# Parameters
# MAXITER=100000
# K0=43
# REGULARISATIONS="prior_on_W L1 prior_on_H L1"
# ITERS=10
# PARALLEL=20
# GROUP=$1
# BASE="/mnt/scratchc/fmlab/smith10/bayesnmf/"
# SCRIPTPATH="${BASE}SignatureAnalyzer-GPU/SignatureAnalyzer-GPU.py"
# INPUTMATRIX="${BASE}sxc/${GROUP}_new_oscn_mixture_tcga_CxS.tsv"
# OUTPUTDIR="${BASE}results/BayesNMF_${GROUP}/"
# python $SCRIPTPATH data $INPUTMATRIX max_iter=$MAXITER output_dir ${OUTPATH} $REGULARISATIONS labeled > ${OUTPATH}/${UUID}_log.txt &
# Code here uses modified version of SignatureAnalyzer-GPU with modification to return data to R session rather than a flat file


## RETICULATE version 1.38
calculateArdNMF <- function(data=NULL,runs=1,max_iter=10000,regularisation_W="L1",regularisation_H="L1",additional.args=NULL,cores=1){
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop(
            "Package \"reticulate\" must be installed to use ARD-NMF provided by SignatureAnalyzer-GPU.",
            call. = FALSE
        )
    }
    # Function to check and setup conda environment - needs more testing
    install_ardnmf()
    # Output directory is not used as object is returned to R but script still
    # requires an argument to be provided so a random 6 character string is used.
    randomOut <- paste0(sample(c(LETTERS,c(0:9)),size = 6,replace = T),collapse = "")
    cmdargs <- list(max_iter=max_iter,output_dir=randomOut,prior_on_W=regularisation_W,prior_on_H=regularisation_H)

    if(!is.null(additional.args)){
        if(is.list(additional.args)){
            cmdargs <- append(cmdargs,additional.args)
        }
    }

    nmf_options <- parse_ard_nmf_args(args = cmdargs)

    if(runs > 1){
        ard_runs <- list()
        for(r in seq_len(runs)){
            cat(paste0("Run ",r," of ",runs,"\n"))
            ard_run <- run_ard_nmf(data = data,nmf_options = nmf_options)
            ard_runs <- append(ard_runs,list(ard_run))
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
    tempout <- main(data = reticulate::r_to_py(tab),arglist = reticulate::r_to_py(nmf_options))
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

choose_ard_nmf <- function(runs = NULL,method="modal",k=NULL,decision="bdiv",specificRun=NULL){
    if(is.null(runs)){
        stop("runs not specified")
    }
    if(!method %in% c("modal","agnostic")){
        stop("unknown selection method")
    }
    if(!decision %in% c("bdiv","obj")){
        stop("unknown decision method")
    }

    if(is.null(specificRun)){
        kvals <- unlist(sapply(runs,FUN = function(x){x$nsigs}))
        if(!is.null(k)){
            if(!k %in% kvals){
                stop("User-specified 'k' number of sigs not present in any run")
            } else {
                message(paste0("selecting user-specified k - ",k))
                nsig <- k
            }
        } else {
            if(method == "modal"){
                nsig <- collapse::fmode(kvals,ties="min")
                message(paste0("selecting modal k - ",nsig))
            } else {
                message("using k-agnostic")
                nsig <- NULL
            }
        }

        if(!is.null(nsig)){
            print(kvals)
            print(nsig)
            print(kvals == nsig)
            subruns <- runs[kvals == nsig]
            message(paste0("runs matching k - ",paste0(names(subruns),collapse = ",")))
        } else {
            subruns <- runs
        }

        best_run <- select_optimal_k(subrun = subruns,decision = decision)

        finalsubruns <- subruns[[best_run]]
    } else {
        if(specificRun %in% names(runs)){
            message(paste0("User-specified run selected -",specificRun))
            finalsubruns <- runs[[specificRun]]
        } else {
            stop("specified run not found in list of runs")
        }
    }
    return(finalsubruns)
}

select_optimal_k <- function(subrun=NULL,decision="bdiv"){
    if(decision == "bdiv"){
        div <- sapply(subrun,FUN = function(x) x$report[nrow(x$report),"b_div"])
    } else if(decision == "obj"){
        div <- sapply(subrun,FUN = function(x) x$report[nrow(x$report),"obj"])
    } else {
        stop("unknown method")
    }
    min_div <- min(div)
    divmin <- which(div %in% min_div)
    best_run_name <- names(div)[divmin]
    message(paste0("best run - ",best_run_name))
    return(best_run_name)
}

# INPUTMATRIX <- "../SignatureAnalyzer-GPU/example_data/POLEMSI_counts_matrix.txt"
# tab <- as.data.frame(read.table("../SignatureAnalyzer-GPU/example_data/POLEMSI_counts_matrix.txt",header = T,sep = "\t",row.names = 1))
# finalrun <- choose_ard_nmf(out10,method = "modal",decision = "bdiv")
#out10 <- calculateArdNMF(data = tab,runs = 20,max_iter = 1000)
# print(finalrun)
