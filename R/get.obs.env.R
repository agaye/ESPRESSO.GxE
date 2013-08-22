#' 
#' @title Generates some observed environmental exposure data
#' @description Adds a set level of error to simulated binary or quantitative data(the true data) to obtain data with a larger variance (the observed data). The level of error is determined by the misclassification rates in binary data and by the set level of variance in the quantitative data.
#' @param env.data A vector of environmental measures that represents the true data
#' @param env.model Model of the exposure: binary=0 ,quantitative-normal=1 or quantitative-uniform=2
#' @param env.prev Prevalence of the environmental determinat
#' @param env.sd Standard deviation
#' @param env.error misclassification rates: 1-sensitivity and 1-specificity
#' @param env.reliability Reliability of the assessment of quantitative exposure
#' @return A dataframe of two coloumns:
#' \code{true.environment} Input data (true data)
#' \code{observed.environment} Observed data
#' @export
#' @author Amadou Gaye

get.obs.env <- function(env.data=NULL,env.model=0,env.prev=0.1,env.sd=1,env.error=c(0.15,0.15),env.reliability=0.8)
{ 

  if(is.null(env.data)){
    cat("\n\n ALERT!\n")
    cat(" No environmental exposure data found.\n")
    cat(" Check the argument 'env.data'\n")
    stop(" End of process!\n\n", call.=FALSE)
   }
   
   true.environment <- env.data
   misclass.rate.1.to.0 <- env.error[1]
   misclass.rate.0.to.1 <- env.error[2]
   numsubs <- length(true.environment)

   if(env.model == 0)
   { # BINARY EXPOSURE (THIS IS EQUIVALENT TO CASE-CONTROL MISCLASSIFICATION 
     # BUT OCCURS AFTER THE SUBJECTS HAVE BEEN SAMPLED INTO THE STUDY)

     # TURN ENV DATA TO "0s" AND "1s" BEFORE USING THE "MISCLASS" FUNCTION
     mean.env <- env.prev
     environ.temp <- env.data + mean.env
         
     # MISCLASSIFY
     environ <- misclassify(environ.temp, misclass.rate.1.to.0, misclass.rate.0.to.1)

     # TURN IT BACK TO ITS INITIAL FORMAT (mean centred)
     observed.environment <- environ - mean.env

   }else{
     if(env.model==1){ 
        var.error <- (env.sd^2/env.reliability)-(env.sd^2)
        # ADD ERROR TO ORIGINAL ENVIRONMENTAL EXPOSURE TO GENERATE OBSERVED DATA
        environ.n <- rnorm(numsubs, true.environment, sqrt(var.error))
     }
     if(env.model==2){ # UNIFORM
        var.error <- (var(env.data)/env.reliability)-(var(env.data))
        env.model.error <- rnorm(numsubs, 0, sqrt(var.error))
        environ.n <- env.data + env.model.error
     }
     observed.environment <- environ.n
   }

   # RETURNED BY THE FUNCTION
   dt <- data.frame(true.environment, observed.environment)
}

