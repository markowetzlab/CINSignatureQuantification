## Code to generate example class objects
# Rerun as required or on class object update

library(CINSignatureQuantification)

data("TCGA_478_Samples_SNP6_GOLD")

t478 <- TCGA_478_Samples_SNP6_GOLD
subsample <- unique(t478$sample)[1:10]
t478 <- t478[t478$sample %in% subsample]

CNobj <- createCNQuant(t478)
CNobj.drews <- calculateFeatures(object = CNobj,method = "drews")
CNobj.drews <- calculateSampleByComponentMatrix(CNobj.drews)
SIGobj.drews <- calculateActivity(CNobj.drews)

CNobj.mac <- calculateFeatures(object = CNobj,method = "mac")
CNobj.mac <- calculateSampleByComponentMatrix(CNobj.mac)
SIGobj.mac <- calculateActivity(CNobj.mac)

usethis::use_data(CNobj.drews,overwrite = TRUE)
usethis::use_data(SIGobj.drews,overwrite = TRUE)
usethis::use_data(CNobj.mac,overwrite = TRUE)
usethis::use_data(SIGobj.mac,overwrite = TRUE)
