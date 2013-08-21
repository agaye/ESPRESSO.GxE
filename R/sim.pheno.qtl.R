#' 
#' @title Simulates continuous outcome data
#' @description Uses the effects data of the determinants to construct a linear predictor(LP). The outcome is normally distributed variable generated with a mean equal to LP and a standard deviation of 1. Some error is then added to the simulated outcome to obtained the observed outcome.
#' @param num.subjects Number of subjects to simulate
#' @param genotype Genotype
#' @param geno.efkt Effects if the genetic variants
#' @param environment Exposure data for environment
#' @param env.efkt Effects of the environmental determiants
#' @param interaction Effect model: main effects=0, Gene-Environment interaction=1, Gene-Gene interaction=2 and Environment-Enviroment interaction=3
#' @param int.efkt Interaction effect
#' @return A dataframe of phenotype
#' @author Amadou Gaye

sim.pheno.qtl <-
function(numsubjects=10000,genotype=NULL,geno.efkt=0.25,environment=NULL,env.efkt=0.5,
         interaction=NULL,int.efkt=1.5)
{  
   # IF GENOTYPE DATA ARE NOT SUPPLIED STOP AND ISSUE AN ALERT
   if(is.null(genotype)){
      cat("\n\n ALERT!\n")
      cat(" No genotype data found.\n")
      cat(" Check the argument 'genotype'\n")
      stop(" End of process!\n\n", call.=FALSE)
   }
   if(is.null(environment)){
      cat("\n\n ALERT!\n")
      cat(" No environmental exposure data found.\n")
      cat(" Check the argument 'environment'\n")
      stop(" End of process!\n\n", call.=FALSE)
   }
   if(is.null(interaction)){
      cat("\n\n ALERT!\n")
      cat(" No interaction data found.\n")
      cat(" This should be the product of 'genotype' by 'environment'\n")
      stop(" End of process!\n\n", call.=FALSE)
   }

   # ALPHA IS EQUAL TO THE MEAN OF THE TRAIT, WHICH IS 0
   num.obs <- numsubjects
   alpha <- 0 

   # GENERATE THE LINEAR PREDICTOR
   lp <- alpha + (geno.efkt*genotype) + (env.efkt*environment) + (int.efkt*interaction)

   # GENERATE THE TRUE PHENOTYPE DATA TO RETURN
   phenotype <- rnorm(num.obs,lp,1) # standardized SD (SD = 1)
   pheno.data <- data.frame(phenotype)
}

