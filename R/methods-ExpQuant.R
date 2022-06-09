setMethod("show", signature=c(object="ExpQuant"),
          definition=function(object){
              cat(class(object),"object\n")
              cat("Experiment name:",object@experimentName,"\n")
              cat("Initialisation:",object@init.date,"\n")
              cat("Last modified:",object@last.modified,"\n")
              cat("Sample count (full):",object@samples.full,"\n")
              cat("Sample count (currrent):",object@samples.current,"\n")
              cat("Genome build:",object@build,"\n")
              cat("Feature method:",object@feature.method,"\n")
          })
