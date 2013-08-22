#' 
#' @title Generates some observed environmental exposure data
#' @description Adds a set level of error to simulated binary or quantitative data(the true data) to obtain data with a larger variance (the observed data).
#' @param phenotype phenotype
#' @param pheno.model Model of the phenotype: binary=0 ,quantitative-normal=1 or quantitative-uniform=2
#' @param pheno.error.1.0 1 to 0 misclassification rate
#' @param pheno.error.0.1 0 to 1 misclassification rate
#' @param pheno.reliability Reliability of the assessment of quantitative phenotype
#' @return A dataframe containing:
#' \code{observed.phenotype} Observed data
#' @export
#' @author Amadou Gaye

get.obs.pheno <-
function (phenotype=NULL, pheno.model, pheno.error.1.0=0.05, pheno.error.0.1=0.05, pheno.reliability=0.9) 
{ 

  if(is.null(phenotype)){
			 cat("\n\n ALERT!\n")
			 cat(" No phenotype data found.\n")
			 cat(" Check the argument 'phenotype'\n")
			 stop(" End of process!\n\n", call.=FALSE)
		} 
		
  # GET THE OBSERVED OUTCOME DATA
  true.phenotype <- phenotype
    
  if(pheno.model==0){ # IF THE OUTCOME IS BINARY
    
    observed.phenotype <- misclassify(true.phenotype, pheno.error.1.0, pheno.error.0.1)
      
  }else{ # IF THE OUTCOME IS CONTINUOUS NORMAL
    
    # USE THE RELIABITLITY OF PHENOTYPE ASSESSMENT TO COMPUTE THE VARIANCE OF MEASURED PHENOTYPES.
    # RELIABITY = (VAR.true.data/VAR.obs.data) AND VAR.obs.data = (VAR.true.data + VAR.measurement)
    # IT FOLLOWS THAT VAR.true.data + VAR.of.estimate = VAR.true.data/RELIABILITY AND THEN:
    # VAR.measurement = (VAR.true.data/RELIABILITY) - VAR.true.data
    var.m <- (1^2/pheno.reliability)-(1^2) # standardized SD (SD = 1)
      
    # GENERATE THE NORMALLY DISTRIBUTED ERROR USING THE ABOVE COMPUTED VARIANCE
    num.obs <- length(phenotype)
    pheno.error <- rnorm(num.obs,0,sqrt(var.m))
      
    # ADD THE ERROR TO ORIGINAL PHENOTYPES TO GENERATE THE OBSERVED PHENOTYPE DATA
    observed.phenotype <- true.phenotype + pheno.error
  } 
  
  # RETURN THE OBSERVED PHENOTYPE DATA AS A DATAFRAME
  df <- data.frame(observed.phenotype)
}
