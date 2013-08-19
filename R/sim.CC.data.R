#' 
#' @title Simulates case and controls
#' @description Generates affected and non-affected subjects until the set sample size is achieved
#' @param num.obs Number of observations to generate per iteration
#' @param numcases Number of cases to simulate
#' @param numcontrols Number of controls to simulate
#' @param allowed.sample.size Maximum number of observations allowed
#' @param disease.prev Prevalence of the binary outcome
#' @param MAF Minor allele frequency
#' @param geno.model Genetic model; 0 for binary and 1 for continuous
#' @param geno.OR Odds ratios of the genetic determinants
#' @param env.model Model of the environmental exposure
#' @param env.prev Prevelance of the environmental determinates
#' @param env.mean Mean under quantitative-normal model
#' @param env.sd Standard deviation under quantitative-normal model
#' @param env.low.lim Lower limit under quantitative-uniform model
#' @param env.up.lim Upper limit under quantitative-uniform model
#' @param env.OR Odds ratios of the environmental determinants
#' @param int.OR   <------------------------------------------------------------------------------------------------
#' @param baseline.OR Baseline odds ratio for subject on 95 percent population centile versus 5 percentile. This parameter reflects the heterogeneity in disease risk arising from determinates that have not been measured or have not been included in the model
#' @return A matrix
#' @author Amadou Gaye

sim.CC.data <-
function(num.obs=20000, numcases=2000, numcontrols=8000, allowed.sample.size=20000000, disease.prev=0.1,
MAF=0.1, geno.model=0, geno.OR=1.5, env.model=0, env.prev=0.1, env.mean=0, env.sd=1, env.low.lim=0, 
env.up.lim=1, env.OR=1.5, int.OR=1.5, baseline.OR=12.36)
{
   # SET UP ZEROED COUNT VECTORS TO DETERMINE WHEN ENOUGH CASES AND CONTROLS HAVE BEEN GENERATED
   complete <- 0
   complete.absolute <- 0
   cases.complete <- 0
   controls.complete <- 0
   block <- 0

   # SET UP A MATRIX TO STORE THE GENERATED DATA
   sim.matrix <- matrix(numeric(0), ncol=6)

   # SET LOOP COUNTER
   numloops <- 0

   # LOOP UNTIL THE SET NUMBER OF CASES AND OR CONTROLS IS ACHIEVED OR THE 
   # THE SET POPULATION SIZE TO SAMPLE FROM IS REACHED
   while(complete==0 && complete.absolute==0)
     {

							# GENERATE THE TRUE GENOTYPE DATA
							geno.data <- sim.geno.data(num.obs, MAF, geno.model)
							allele.A <- geno.data$allele.A
							allele.B <- geno.data$allele.B
							genotype <- geno.data$genotype
							
							# GENERATE THE TRUE ENVIRONMEANTAL EXPOSURE DATA
       environment <- sim.env.data(num.obs, env.model, env.prev, env.mean, env.sd, env.low.lim, env.up.lim)
       
							# GENERATE THE TRUE INTERACTION DATA
       interaction <- genotype*environment       

					  # GENERATE SUBJECT EFFECT DATA THAT REFLECTS BASELINE RISK: 
					  # NORMALLY DISTRIBUTED RANDOM EFFECT VECTOR WITH APPROPRIATE 
					  # VARIANCE ON SCALE OF LOG-ODDS
					  subject.effect.data <- sim.subject.data(num.obs, baseline.OR)
					  
       # GENERATE THE TRUE OUTCOME DATA
       pheno.data <- sim.pheno.bin(num.obs, disease.prev, genotype, environment, interaction, subject.effect.data, 
                                   geno.OR, env.OR, int.OR)
       phenotype <- pheno.data$phenotype

       # STORE THE TRUE OUTCOME, GENETIC AND ENVIRONMENT AND ALLELE DATA IN AN OUTPUT MATRIX 
       # WHERE EACH ROW HOLDS THE RECORDS OF ONE INDIVUDAL
       sim.matrix.temp <- cbind(phenotype,genotype,allele.A,allele.B,environment,interaction)

       # UPDATE THE MATRIX THAT HOLDS ALL THE DATA GENERATED SO FAR, AFTER EACH LOOP
       sim.matrix <- rbind(sim.matrix, sim.matrix.temp)

       # SELECT OUT CASES
       sim.matrix.cases <- sim.matrix[phenotype==1,]

       # SELECT OUT CONTROLS
       sim.matrix.controls <- sim.matrix[phenotype==0,]

       # COUNT THE NUMBER OF CASES AND CONTROLS THAT HAS BEEN GENERATED
       cases.simulated <- dim(sim.matrix.cases)[1]
       controls.simulated <- dim(sim.matrix.controls)[1]

       # TEST IF THERE ARE AT LEAST ENOUGH CASES ALREADY SIMULATED
       # IF THERE ARE, DEFINE THE CASE ELEMENT OF THE DATA MATRIX
       if(cases.simulated >= numcases)
       {
         sim.matrix.cases <- sim.matrix.cases[1:numcases,]
         cases.complete <- 1
       }

       # TEST IF THERE ARE AT LEAST ENOUGH CONTROLS ALREADY SIMULATED
       # IF THERE ARE, DEFINE THE CONTROL ELEMENT OF THE DATA MATRIX
       if(controls.simulated>=numcontrols)
       {
         sim.matrix.controls <- sim.matrix.controls[1:numcontrols,]
         controls.complete <- 1
       }

       # HAVE WE NOW GENERATED THE SET NUMBER OF CASES AND CONTROLS?
       complete <- cases.complete*controls.complete		

       # HAVE WE EXCEEDED THE TOTAL SAMPLE SIZE ALLOWED?
       complete.absolute <- (((block+1)*num.obs)>=allowed.sample.size)
       if(complete.absolute==1) {sample.size.excess <- 1}else{sample.size.excess <- 0}

        # INCREMENT LOOP COUNTER
        numloops <- numloops + 1
    }

   # STACK FINAL DATA MATRIX WITH CASES FIRST
   sim.matrix <- rbind(sim.matrix.cases,sim.matrix.controls)
   totalnumrows <- dim(sim.matrix)[1]
   sim.matrix <- cbind(1:totalnumrows, sim.matrix)

   # NAME THE COLUMNS OF THE MATRIX AND RETURN IT AS A DATAFRAMEDATAFRAME
   colnames(sim.matrix) <- c("id", "phenotype", "genotype", "allele.A", "allele.B", "environment", "interaction")
   mm <- list(data=data.frame(sim.matrix), allowed.sample.size.exceeded=sample.size.excess)
}

