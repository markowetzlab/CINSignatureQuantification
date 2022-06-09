## Testing ground for package and functions

# - CHECK E/W/N
#   - Non-standard license
# - pkgdown
# - clinPredictionPlatinum-methods uses normalised but not threshold adjusted

# Load library
library(CINSignatureQuantification)

# Load test data
dfTest = readRDS("inst/TCGA_478_Samples_SNP6_GOLD.rds")

## Pipeline method
sigAct478.drews = quantifyCNSignatures(dfTest,experimentName = "478TCGAPCAWG",method = "drews",cores = 6)
sigAct478.mac = quantifyCNSignatures(dfTest,experimentName = "478TCGAPCAWG",method = "mac",cores = 6)

## Individual functions
# Convert to CNQuant object
myData = createCNQuant(data = dfTest)

## Feature extraction (includes smoothing and preparing data)
myData.drews = calculateFeatures(myData, method="drews",cores = 1)
myData.mac = calculateFeatures(myData, method="mac",cores = 1)

## Get sum-of-posterior matrix
myData.drews = calculateSampleByComponentMatrix(myData.drews)
myData.mac = calculateSampleByComponentMatrix(myData.mac)

## Get activities
myData.drews = calculateActivity(myData.drews)
myData.mac = calculateActivity(myData.mac)

## Test clinical classifier (CX3/CX2 and De-novo for two self chosen signatures)
vPredPlat = clinPredictionPlatinum(sigAct478.drews)
vPredCX8CX9 = clinPredictionDenovo(sigAct478.drews, sampTrain = sample(getSamples(sigAct478.drews), 50), sigsTrain = c("CX9", "CX8"))

## Additional functions
# Show and subsetting
sigAct478.drews
sigAct478.drews[1:10]
sigAct478.drews[getSamples(sigAct478.drews)[1:50]]

# Sample feature/clinical information
getSamplefeatures(sigAct478.drews)
load("inst/test.sample.features.rda")
sigAct478.drews = addsampleFeatures(object = sigAct478.drews,sample.data = test.sample.features)
getSamplefeatures(sigAct478.drews)

# plots
plotSampleByComponent(sigAct478.drews)
plotSegments(sigAct478.drews,sample = 1,cn.max = 8)
plotActivities(object = sigAct478.drews)

# misc
getSampleByComponent(sigAct478.drews)
getExperiment(sigAct478.drews)
getSamples(sigAct478.drews)
