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

INPUTMATRIX <- "../SignatureAnalyzer-GPU/example_data/POLEMSI_counts_matrix.txt"
tab <- as.data.frame(read.table("../SignatureAnalyzer-GPU/example_data/POLEMSI_counts_matrix.txt",header = T,sep = "\t",row.names = 1))

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

    finalRun <- choose_ard_nmf(ard_runs)
    print(finalRun)
    return(ard_runs)
}

choose_ard_nmf <- function(runs = NULL){
    modalSigs <- collapse::fmode(unlist(lapply(runs,FUN = function(x){ x$nsigs})))
    return(modalSigs)
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
    y$report <- x[[2]]
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

out10 <- calculateArdNMF(data = tab,runs = 10,max_iter = 1000)
