#'
#' @title Generates exposure data with some error
#' @description Uses functions make.obs.geno and make.obs.env to generate effect data with a set level of error
#' @param true.data Input table of simulated data considered as true data
#' @param geno.error Misclassification rates in genetic assessment: 1-sensitivity and 1-specificity
#' @param geno.model Genetic model; 0 for binary and 1 for continuous
#' @param MAF Minor allele frequency
#' @param env.error Misclassification rates in environmental exposures assessment: 1-sensitivity and 1-specificity
#' @param env.model Model of the exposure: binary=0, quantitative-normal=1 or quantitative-uniform=2
#' @param env.prev Prevalence of environmental exposure
#' @param env.sd Standard Deviation
#' @param env.reliability Reliability of the assessment of quantitative exposure
#' @param pheno.model Model of the phenotype: binary=0 ,quantitative-normal=1 or quantitative-uniform=2
#' @param pheno.error Misclassification rate: 1 to 0, 0 to 1
#' @param pheno.reliability Reliability of the assessment of quantitative phenotype
#' @return A matrix
#' @export
#' @author Amadou Gaye

get.observed.data <-
function(true.data=NULL,geno.error=c(0.05,0.05),geno.model=0,MAF=0.1,
         env.error=c(0.15,0.15),env.model=0,env.prev=0.1,env.sd=1,env.reliability=0.8,
         pheno.model=0,pheno.error=c(0.1,0.1),pheno.reliability=0.9)
{
    if(is.null(true.data)){
      cat("\n\n ALERT!\n")
      cat(" No data found.\n")
      cat(" Check the argument 'true.data'\n")
      stop(" End of process!\n\n", call.=FALSE)
    }
		 	
		 	sim.df <- true.data      

    # GET THE OBSERVED GENOTYPES
    geno.error.1.0 <- geno.error[1]
    geno.error.0.1 <- geno.error[2]
    pheno.error.1.0 <- pheno.error[1]
    pheno.error.0.1 <- pheno.error[2]
    true.genotype <- sim.df$genotype
    obs.genotype <- get.obs.geno(sim.df$allele.A,sim.df$allele.B, geno.model, geno.error.1.0,geno.error.0.1,MAF)
    
    # GET THE OBSERVED ENVIRONMENTAL EXPOSURE DATA
    true.environment <- sim.df$environment
    obs.environment <- get.obs.env(true.environment,env.model,env.prev,env.sd,env.error,env.reliability)
    
    # GET THE OBSERVED INTERACTION DATA
    obs.interaction <- obs.genotype$observed.genotype*obs.environment$observed.environment

    # GET THE OBSERVED OUTCOME DATA    
    true.phenotype <- sim.df$phenotype
    obs.phenotype <- get.obs.pheno(true.phenotype,pheno.model,pheno.error.1.0,pheno.error.0.1,pheno.reliability)

    # REPLACE THE TRUE DATA BY THE NOW GENERATED OBSERVED GENOTYPES
    # IN THE INITIAL MATRIX THAT HELD THE TRUE DATA
    sim.df$genotype <- obs.genotype$observed.genotype
    sim.df$allele.A <- obs.genotype$observed.allele.A
    sim.df$allele.B <- obs.genotype$observed.allele.B
    sim.df$environment <- obs.environment$observed.environment
    sim.df$phenotype <- obs.phenotype$observed.phenotype
    
    # RETURN THE MATRIX WHICH NOW CONTAINS ONLY THE OBSERVED DATA TO ANALYSE BY GLM
    colnames(sim.df) <- c("id", "phenotype", "genotype", "allele.A", "allele.B", "environment", "interaction")
    return(sim.df)
}

