#'
#' @title Simulates subjects for continuous outcome
#' @description Generates the specified number of subjects
#' @param num.subjects Number of subjects to simulate
#' @param MAF Minor allele frequencies of the two genetic variants
#' @param geno.model Genetic model; 0 for binary and 1 for continuous
#' @param geno.efkt Effects of the genetic variants
#' @param env.model Model of the environmental exposure
#' @param env.efkt Effects of the environment determinats
#' @param env.prev Prevalences of the environmental determinats
#' @param env.mean Mean under quantitative-normal model
#' @param env.sd Standard deviation under quantitative-normal model
#' @param env.low.lim lower limit under quantitative-uniform model
#' @param env.up.lim upper limit under quantitative-uniform model
#' @param int.efkt Interaction effect
#' @return A matrix
#' @author Amadou Gaye

sim.QTL.data <-
function(numsubjects=500,MAF=0.1,geno.model=0,geno.efkt=0.25,env.model=0,env.efkt=0.5,
         env.prev=0.1,env.mean=0,env.sd=1,env.low.lim=0,env.up.lim=1,int.efkt=1.5)
{
   num.obs <- numsubjects
	   
   # GENERATE THE TRUE GENOTYPE DATA
   geno.data <- sim.geno.data(num.obs, MAF, geno.model)
   allele.A <- geno.data$allele.A
   allele.B <- geno.data$allele.B
   genotype <- geno.data$genotype
			
   # GENERATE THE TRUE ENVIRONMENTAL EXPOSURE DATA			
   environment <- sim.env.data(num.obs,env.model,env.prev,env.mean,env.sd,env.low.lim,env.up.lim)
          
   # GENERATE THE TRUE INTERACTION DATA           
   interaction <- genotype*environment
           
   # GENERATE OUTCOME DATA
   pheno.data <- sim.pheno.qtl(num.obs,genotype,geno.efkt,environment,env.efkt,interaction,int.efkt)
   phenotype <- pheno.data$phenotype

   # STORE THE GENERATED TRUE DATA INTO AN OUTPUT MATRIX 
   sim.matrix <- cbind(phenotype,genotype,allele.A,allele.B,environment,interaction)

   # ADD IDs (JUST A ROW COUNT)
   totalnumrows <- dim(sim.matrix)[1]
   sim.matrix <- cbind(1:totalnumrows, sim.matrix)

   # ADD COLUMN NAMES AND RETURN A DATAFRAME
   colnames(sim.matrix) <- c("id","phenotype","genotype","allele.A","allele.B","environment","interaction")
   mm <- data.frame(sim.matrix)
}

