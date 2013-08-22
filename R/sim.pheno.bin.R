#'
#' @title Generates phenotype status
#' @description Generates affected and non-affected subjects
#' @param num.obs Number of observations to generate per iteration
#' @param disease.prev Prevalence of the binary outcome
#' @param genotype Exposure data for genetic determinate
#' @param environment Exposure data for environment
#' @param interaction Effect model: main effects=0, Gene-Environment interaction=1, Gene-Gene interaction=2 and Environment-Enviroment interaction=3
#' @param subject.effect.data Subject effect data, reflects the heterogenity in baseline disease risk
#' @param geno.OR Odds ratios of the two genetic determinates
#' @param env.OR Odds ratios of the two environments
#' @param int.OR int.OR Odds ration of the interaction
#' @return A dataframe of phenotype
#' @author Amadou Gaye

sim.pheno.bin <-
function(num.obs=20000, disease.prev=0.1, genotype=NULL, environment=NULL, interaction=NULL, subject.effect.data=NULL, 
         geno.OR=1.5, env.OR=1.5, int.OR=1.5)
{ 
   # IF GENOTYPE AND SUBJECT EFFECT DATA ARE NOT SUPPLIED STOP AND ISSUE AN ALERT
   if(is.null(genotype)){
			  cat("\n\n ALERT!\n")
			  cat(" No genotype data found.\n")
			  cat(" Check the argument 'genotype'\n")
			  stop(" End of process!\n\n", call.=FALSE)
			}
			if(is.null(environment)){
			  cat("\n\n ALERT!\n")
			  cat(" No environment data found.\n")
			  cat(" Check the argument 'environment'\n")
			  stop(" End of process!\n\n", call.=FALSE)
			}
		 if(is.null(interaction)){
			  cat("\n\n ALERT!\n")
			  cat(" No interaction data found.\n")
			  cat(" This should be the product 'genotype' by 'environment'\n")
			  stop(" End of process!\n\n", call.=FALSE)
			}
			if(is.null(subject.effect.data)){
			  cat("\n\n ALERT!\n")
			  cat(" No baseline effect data found.\n")
			  cat(" Check the argument 'subject.effect.data'\n")
			  stop(" End of process!\n\n", call.=FALSE)
			}
 
   # GET THE ALPHA AND BETA VALUES
   alpha <- log(disease.prev/(1-disease.prev))
   geno.beta <-	log(geno.OR)
   env.beta <-	log(env.OR)
   int.beta <- log(int.OR)

   # GENERATE THE LINEAR PREDICTOR
   lp <- alpha + (geno.beta*genotype) + (env.beta*environment) + (int.beta*interaction) + subject.effect.data
   # GET THE 'mu' THE PROBABILITY OF DISEASE THROUGH LOGISTIC TRANSFORMATION
   mu <- exp(lp)/(1 + exp(lp))
   
   # GENERATE THE PHENOTYPE DATA AND RETURN IT AS A DATAFRAME
   phenotype <- rbinom(num.obs,1,mu)
   pheno.data <- data.frame(phenotype)
}

