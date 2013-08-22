#'
#'@title Adds some error to genotype data
#'@description Simulates errors and adds it to the true data to obtain observed data. The alleles simulated by the function sim.geno.data are randomly misclassified and used to form new genotypes that represent the observed genotypes.
#'@param allele.A Allele A
#'@param allele.B Allele B
#'@param geno.model Genetic model; 0 for binary and 1 for continuous
#'@param MAF Minor Allele Frequency
#'@param geno.error.0.1 0 to 1 misclassififcation rate
#'@param geno.error.1.0 1 to 0 misclassififcation rate
#'@return A dataframe containing:
#'\code{observed.genotype} Observed genotypes
#'\code{observed.allele.A} Observed A alleles
#'\code{observed.allele.B} Observed B alleles
#'@export
#'@author Amadou Gaye

get.obs.geno <-
function (allele.A=NULL, allele.B=NULL, geno.model=0, MAF=0.1, geno.error.0.1=0.05, geno.error.1.0=0.05) 
{
   # IF ALLELE DATA ARE NOT SUPPLIED STOP AND ISSUE AN ALERT
   if(is.null(allele.A)){
			  cat("\n\n ALERT!\n")
			  cat(" No allele data found for the first allele.\n")
			  cat(" Check the argument 'allele.A'\n")
			  stop(" End of process!\n\n", call.=FALSE)
			}
			if(is.null(allele.B)){
			  cat("\n\n ALERT!\n")
			  cat(" No allele data found for the second allele.\n")
			  cat(" Check the argument 'allele.B'\n")
			  stop(" End of process!\n\n", call.=FALSE)
			}
			
   mean.add <- (2*MAF*(1-MAF)+2*(MAF^2))
   mean.bin <- (2*MAF*(1-MAF)+(MAF^2))
   true.allele.A <- allele.A
   true.allele.B <- allele.B

   observed.allele.A <- misclassify(true.allele.A, geno.error.1.0, geno.error.0.1)
   observed.allele.B <- misclassify(true.allele.B, geno.error.1.0, geno.error.0.1)
   observed.genotype <- observed.allele.A+observed.allele.B

   if(geno.model==0){
      genotyp.U <- observed.genotype > 0
      observed.genotype <- genotyp.U - mean.bin
   }
   if(geno.model==1){
      genotyp.U <- observed.genotype
      observed.genotype <- genotyp.U - mean.add
   }

   df <- data.frame(observed.genotype, observed.allele.A, observed.allele.B)
}

