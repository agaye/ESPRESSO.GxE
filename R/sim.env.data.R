#' 
#' @title Simulates cases and controls
#' @description Generates data for a binary, quantitative-normal, or quantitative-uniform environmental determinate
#' @param num.obs Number of observations to simulate
#' @param env.model Model of the exposure: binary=0, quantitative-normal=1, quantitative-uniform=2
#' @param env.prev Prevalence of the environmental exposure
#' @param env.mean Mean under quantitative-normal model
#' @param env.sd Standard deviation under quantitative-normal model
#' @param env.low.lim Lower limit under quantitative-uniform model
#' @param env.up.lim Upper limit under quantitative-uniform model
#' @return A vector of continuous or binary values
#' @author Amadou Gaye

sim.env.data <-
function (num.obs=20000,env.model=0,env.prev=0.1,env.mean=0,env.sd=1,env.low.lim=0,env.up.lim=1) 

{ 	
   # CREATE THE FIRST ENVIRONMENTAL COVARIATE 
   if(env.model==0){   # BINARY DISTRIBUTION
      env.U <- rbinom(num.obs, 1, env.prev)
      mean.e <- env.prev
      env.U <- env.U-mean.e
   }
   if(env.model==1){   # NORMAL DISTRIBUTION
      env.U <- rnorm(num.obs, env.mean, env.sd)
      env.U <- env.U-mean(env.U)     # mean centering
   }
   if(env.model==2){  # UNIFORM DISTRIBUTION
      if(env.low.lim >= env.up.lim)
      {
        cat("\n\nALERT!\n Uniform Distribution: The upper limit must be greater than the lower limit\n\n")
      }
      env.U <- runif(num.obs, env.low.lim, env.up.lim)
      env.U <- env.U-mean(env.U)     # mean centering
   }

   # return a vector 
   environment <- env.U
}

