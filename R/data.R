#' gap_hg19
#'
#' Chromosomal banding and position of genomic features in genome build hg19
#'
#' @docType data
#' @keywords datasets
#' @name gap_hg19
#' @usage data(gap_hg19)
#' @format A data frame with 457 rows and 9 variables
NULL

#' hg19.chrom.sizes
#'
#' Chromosomal lengths for genome build hg19
#'
#' @docType data
#' @keywords datasets
#' @name hg19.chrom.sizes
#' @usage data(hg19.chrom.sizes)
#' @format A data frame with 24 rows and 2 variables
NULL

#' gap_hg38
#'
#' Chromosomal banding and position of genomic features in genome build hg38
#'
#' @docType data
#' @keywords datasets
#' @name gap_hg38
#' @usage data(gap_hg38)
#' @format A data frame with 827 rows and 9 variables
NULL

#' hg38.chrom.sizes
#'
#' Chromosomal lengths for genome build hg38
#'
#' @docType data
#' @keywords datasets
#' @name hg38.chrom.sizes
#' @usage data(hg38.chrom.sizes)
#' @format A data frame with 24 rows and 2 variables
NULL

#' Drews2022_CX3CX2_Clinical_classifier
#'
#' List of mean and scaling factors for signatures 2 and 3 for the prediction of
#' platinum status used in Drews 2022
#'
#' @docType data
#' @keywords datasets
#' @name Drews2022_CX3CX2_Clinical_classifier
#' @usage data(Drews2022_CX3CX2_Clinical_classifier)
#' @format A list of 2 numeric vectors of length 2
NULL

#' Drews2022_TCGA_Mixture_Models
#'
#' List of mixture model mean, standard deviation, and weight for each copy number
#' feature used in Drews 2022
#'
#' @docType data
#' @keywords datasets
#' @name Drews2022_TCGA_Mixture_Models
#' @usage data(Drews2022_TCGA_Mixture_Models)
#' @format A list of 5 data.frames
NULL

#' Drews2022_TCGA_Scaling_Variables
#'
#' List of mean and scaling factors for signatures used in Drews 2022
#'
#' @docType data
#' @keywords datasets
#' @name Drews2022_TCGA_Scaling_Variables
#' @usage data(Drews2022_TCGA_Scaling_Variables)
#' @format A list of 2 numeric vectors of length 17
NULL

#' Drews2022_TCGA_Signature_Thresholds
#'
#' Numeric vector containing signature-specific thresholds used in activity
#' calculations used in Drews 2022
#'
#' @docType data
#' @keywords datasets
#' @name Drews2022_TCGA_Signature_Thresholds
#' @usage data(Drews2022_TCGA_Signature_Thresholds)
#' @format A numeric vector of length 17
NULL

#' Drews2022_TCGA_Signatures
#'
#' Signature-by-component matrix for 17 derived signatures as used in Drews 2022
#'
#' @docType data
#' @keywords datasets
#' @name Drews2022_TCGA_Signatures
#' @usage data(Drews2022_TCGA_Signatures)
#' @format A 17 by 43 numeric matrix
NULL

#' Macintyre2018_OV_Mixture_Models
#'
#' List of mixture model mean, standard deviation, and weight for each copy number
#' feature used in Macintyre 2018
#'
#' @docType data
#' @keywords datasets
#' @name Macintyre2018_OV_Mixture_Models
#' @usage data(Macintyre2018_OV_Mixture_Models)
#' @format A list of 6 data.frames
NULL

#' Macintyre2018_OV_Signatures
#'
#' Signature-by-component matrix for 7 derived signatures as used in Macintyre 2018
#'
#' @docType data
#' @keywords datasets
#' @name Macintyre2018_OV_Signatures
#' @usage data(Macintyre2018_OV_Signatures)
#' @format A 7 by 36 numeric matrix
NULL

#' Macintyre2018_OV_Signatures_normalised
#'
#' Signature-by-component matrix for 7 derived signatures as used in Macintyre 2018
#' where column sums are normalised to 1.
#'
#' @docType data
#' @keywords datasets
#' @name Macintyre2018_OV_Signatures_normalised
#' @usage data(Macintyre2018_OV_Signatures_normalised)
#' @format A 7 by 36 numeric matrix
NULL

#' cancerSpecificSignatures
#'
#' List of 20 named character vectors describing the signatures associated with each
#' cancer type (as described in Drews et al. 2022 signature discovery). Vectors include
#' 10 pan-cancer signatures plus any additional signatures specific to the selected cancer type.
#'
#' Cancer types are limited to "BRCA", "BLCA", "OV", "LUAD", "GBM", "STAD", "SKCM", "UCEC", "HNSC",
#' "CESC", "LUSC", "SARC", "COAD", "KIRC", "LGG", "READ", "PRAD", "ESCA", "LIHC", & "TGCT"
#'
#' @docType data
#' @keywords datasets
#' @name cancerSpecificSignatures
#' @usage data(cancerSpecificSignatures)
#' @format A length 20 named list of character vectors
NULL

#' TCGA_478_Samples_SNP6_GOLD
#'
#' Example data set of 478 pan-cancer samples with detectable CIN in the correct
#' input format for generating copy number signatures. Data was generated during
#' the development and submission of "A pan-cancer compendium of CIN" drews et
#' al. 2022.
#'
#' @docType data
#' @keywords datasets
#' @name TCGA_478_Samples_SNP6_GOLD
#' @usage data(TCGA_478_Samples_SNP6_GOLD)
#' @format A data.frame with 5 columns and 65098 rows.
NULL

#' test.sample.features
#'
#' Example data set of additional example sample features to add to the existing
#' sample feature data of the 478 pan-cancer samples with detectable CIN. Data
#' in this dataset are dummy variables only for function testing and examples
#' and have no meaning.
#'
#' @docType data
#' @keywords datasets
#' @name test.sample.features
#' @usage data(test.sample.features)
#' @format A data.frame with 4 columns and 478 rows.
NULL
