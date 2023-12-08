# addsampleFeatures

data("CNobj.drews")
data("SIGobj.drews")
data("test.sample.features")

cnTrue <- CNobj.drews
cnTrue@samplefeatData <- merge.data.frame(cnTrue@samplefeatData,
                                          test.sample.features,
                                          by.x = "row.names",
                                          by.y = "sample",all.x = T)
rownames(cnTrue@samplefeatData) <- cnTrue@samplefeatData$Row.names
cnTrue@samplefeatData <- cnTrue@samplefeatData[,-1]

sigTrue <- SIGobj.drews
sigTrue@samplefeatData <- merge.data.frame(sigTrue@samplefeatData,
                                          test.sample.features,
                                          by.x = "row.names",
                                          by.y = "sample",all.x = T)
rownames(sigTrue@samplefeatData) <- sigTrue@samplefeatData$Row.names
sigTrue@samplefeatData <- sigTrue@samplefeatData[,-1]

test_that("test null", {
    expect_error(addsampleFeatures())
})

test_that("test class input", {
    expect_error(addsampleFeatures(object=test.sample.features,
                                   sample.data=test.sample.features))
})

test_that("test missing sample.data", {
    expect_error(addsampleFeatures(object=CNobj.drews))
})

test_that("test sample.data type", {
    expect_error(addsampleFeatures(object=CNobj.drews,
                                   sample.data=test.sample.features$age))
})

test_that("test CNobj correct", {
    expect_equal(addsampleFeatures(object=CNobj.drews,
                                   sample.data=test.sample.features)@samplefeatData,
                 cnTrue@samplefeatData)
})

test_that("test SIGobj correct", {
    expect_equal(addsampleFeatures(object=SIGobj.drews,
                                   sample.data=test.sample.features)@samplefeatData,
                 sigTrue@samplefeatData)
})

test.sample.features$sample <- paste0("sample_",seq.int(1,nrow(test.sample.features),1))
test_that("test no matches", {
    expect_error(addsampleFeatures(object=CNobj.drews,
                                   sample.data=test.sample.features))
})

colnames(test.sample.features) <- c("samps","type","age","location")
test_that("test no matches", {
    expect_error(addsampleFeatures(object=CNobj.drews,
                                   sample.data=test.sample.features))
})
