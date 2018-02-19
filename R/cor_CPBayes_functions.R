 ##=================================********** Correlated version of CPBayes ***********==============================##
 ## This function is the main MCMC function for implementing CPBayes in case of correlated summary statistics.
 ## Correlated summary statistics arise for case-control statistics with overlapping subjects or cohort data. 
 ## It calls the function: initiate_MCMC() from 'common_MCMC_functions.R' to initialize the parameters in the MCMC.
 ## It calls: correlated_beta_update() from 'cor_MCMC_functions.R'; & Z_update(), q_update(), de_update1(), de_update0() 
 ## from 'common_MCMC_functions.R' to update different parameters in MCMC.
 ## It calls select_subset(), overall_pleio_measure() from 'summary_functions.R' to summarize the MCMC data.
 ##=================================***************************************************===============================## 

 ##=================================**************Argument explainations***************===============================##
 ## variantName: name of the genetic variant. It must be a character vector of length of 1.
 ## traitNames:  name of phenotypes. It must be a character vector containing the phenotypes names. 
 ## X:           beta hat vector. 
 ## S:           covariance matrix of X (beta hat) for the phenotypes.
 ## updateDE:    logical. Indicates whether to update the parameter in MCMC or not.  
 ## RP:          total number of replications in the MCMC. 
 ## burn.in:     burn in period - after which the MCMC sample should be selected.
 ## Note that, traitNames, X, s.e., and S, all will have the same order in phenotypes. 
 ##=================================***************************************************===============================## 

 ##library("MASS")                                               ## MASS package is required

 CPBayes_cor = function(variantName, traitNames, X, s.e., corln, updateDE, MinSlabVar, MaxSlabVar, RP, burn.in)
 {
     set.seed(10)
     K = length(X)
     PPAj_thr = 0.20                                           ## PPAj threshold
     
     ## If the covariance matrix is not positive definite, the diagonal elements are incremented to make it PD 
     S <- diag(s.e.)%*%corln%*%diag(s.e.)
     epsilon = 10^(-5)
     increment = rep(epsilon,K)
     while(det(S) <= 0)
     {   diag(S) = diag(S)+increment
     }
     S.inv = solve(S)                                          ## inverse of the cov matrix to be passed in MCMC function

     ## set the initial choices of parameters     
     tau <- 0.01                                               ## spike sd
     CentralSlabVar <- (MinSlabVar+MaxSlabVar)/2
     nonNullVar <- CentralSlabVar                              ## central value of slab variance 
     de <- sqrt(tau^2/nonNullVar)                              ## 1/de - ratio of spike and slab variances

     ## specfications for updating 'de' parameter
     min_var <- MinSlabVar                                     ## minimum value of spike variance
     max_var <- MaxSlabVar                                     ## maximum value of slab variance
     max_de <- tau/sqrt(min_var)                               ## maximum choice of 'de' parameter
     min_de <- tau/sqrt(max_var)                               ## minimum choice of 'de'
     shape1_de <- 1                                            ## choice of shape1 parameter of Beta(shape1,1) prior of updating 'de'

     ## an informed initialization 
     initiate <- initiate_MCMC( K, X, s.e. )
     beta <- initiate$beta
     Z <- initiate$Z
     q <- initiate$q
     nA <- initiate$K1_FDR                                     ## number of associated traits

     ## specify the shape parameters of 'q' prior
     shape2 <- 1                                               ## shape2 parameter for Beta prior of q
     LB <- 0.1; UB <- 0.5;
     qm <- nA/K
     if(qm < LB) qm <- LB
     if(qm > UB) qm <- UB

     #qm <- 0.25
     shape1 <- (qm/(1-qm)) * shape2

     thinning <- 1                                             ## thinning period for MCMC


     mcmc.samplesize <- (RP-burn.in) %/% thinning; Z.data <- matrix(0,mcmc.samplesize,K); row <- 0 ;
     sim.beta <- matrix(0,mcmc.samplesize,K); sample_probZ_zero <- matrix(0,mcmc.samplesize,K);

     for(rp in  1:RP)
     {
         ## Update the main beta parameters                                              
         beta = correlated_beta_update(K, tau, de, Z, X, S.inv)

         ## Update Z
         res_Z <- Z_update(K, q, tau, de, beta)
         Z <- res_Z$Z
         probZ_zero <- res_Z$prob
         ## Update q
         q = q_update(K, Z, shape1, shape2)

         ## Update de
         if(updateDE == TRUE){
             if(sum(Z) > 0) 
                de <- de_update1(min_de, max_de, shape1_de, beta, Z, tau) 
             else de <- de_update0(min_de, max_de, shape1_de)
         }

         if(rp > burn.in && rp%%thinning == 0)
         {   row = row + 1 ;
             Z.data[row,] = Z ;
             sim.beta[row,] = beta ;
             sample_probZ_zero[row,] = probZ_zero ;
         }

     } # close rp      

     ##......................................... Compute the summary .................................................##
     ## selection of subset 
     cor.subset <- select_subset( K, Z.data, mcmc.samplesize )
     selected_traits <- NULL
     if( length(cor.subset) > 0 ) selected_traits <- traitNames[cor.subset]

     ## extracting traits having PPAj > PPAj_thr         
     asso.pr = colSums(Z.data)/mcmc.samplesize
     which_traits = which(asso.pr > PPAj_thr)
     
     imp_PPAj = 0; imp_traits = 0; important_phenos = NULL
     
     if(length(which_traits) > 0){
     	imp_PPAj = asso.pr[which_traits]
     	imp_traits = traitNames[which_traits]
     	important_phenos = data.frame( traits = imp_traits, PPAj = imp_PPAj, stringsAsFactors = FALSE)
     } 
     
     ## calculate the Bayes factor and PPNA
     pleio_evidence <- overall_pleio_measure( K, shape1, shape2, sample_probZ_zero )
     log10_BF_cor <- pleio_evidence$log10_BF
     PPNA.cor <- pleio_evidence$PPNA    

     ## return the outputs
     data = list( variantName = variantName, log10_BF = log10_BF_cor, locFDR = PPNA.cor, subset = selected_traits, 
                  important_traits = important_phenos, auxi_data = list( traitNames = traitNames, K = K, 
                  mcmc.samplesize = mcmc.samplesize, PPAj = asso.pr, Z.data = Z.data, sim.beta = sim.beta, betahat = X, se = s.e. ) )
 }
 
 

 ##=================================***********Combined strategy of CPBayes************===============================##
 ## This function combines the uncorrelated and correlated versions of CPBayes to propose a combined strategy.
 ## It first runs the correlated version of CPBayes. 
 ## Then check whether the phenotypes having smallest univariate p-values are selected. 
 ## If the checking answers 'negative', we run the uncorrelated version and accept the results obtained.
 ## So, it first calls CPBayes_cor(), then if necessary, it calls CPBayes_uncor()  
 ## This is the primary function for performing correlated version of CPBayes
 ##=================================***************************************************===============================## 

 ##=================================**************Argument explainations***************===============================##
 ## variantName: name of the genetic variant. It must be a character vector of length of 1
 ## traitNames:  name of phenotypes. It must be a character vector containing the phenotypes names 
 ## X:           beta hat vector 
 ## s.e.:        standard error vector 
 ## corln:       correlation matrix of X (beta hat) for the phenotypes
 ## updateDE:    logical. Indicates whether to update the parameter in MCMC or not  
 ## RP:          total number of replications in the MCMC 
 ## burn.in:     burn in period, after which the MCMC sample should be selected
 ## Note that, traitNames, X, s.e., and corln, all will have the same order in phenotypes 
 ##=================================***************************************************===============================## 

 ## combined strategy using both of the correlated and uncorrelated versions of the Bayes meta analysis for using in correlated case
 
 combined_CPBayes = function(variantName, traitNames, X, s.e., corln, updateDE, MinSlabVar, MaxSlabVar, RP, burn.in)
 {
     ptm1 <- proc.time()
     K = length(X)
     
     ## First, run correlated version of CPBayes
     res_cor = CPBayes_cor(variantName, traitNames, X, s.e., corln, updateDE, MinSlabVar, MaxSlabVar, RP, burn.in)           

     selected_traits = res_cor$subset
     cor.subset <- match(selected_traits, traitNames)
     ## compute the length of the subset of traits selected by CPBayes
     K1_cor = length(cor.subset)     

     if(K1_cor > 1) cor.subset = sort(cor.subset)

     indi_uncor = 0                                            ## whether the uncorrelated version of CPBayes was chosen for final analysis

     if(K1_cor > 0){
         pv = pchisq( (X/s.e.)^2, df=1, lower.tail=F )
         sort.pv = sort(pv, index.return = TRUE)
         sorted = sort.pv$x 
         sorted.index = sort.pv$ix 
 
         pv.subset = sorted.index[1:K1_cor]
         pv.subset = sort(pv.subset)
         ## checking whether correlated CPBayes selected traits having smalled p-values
         diff = pv.subset-cor.subset
         Sum = diff%*%diff
         ## implementing uncor CPBayes if the checking is not satisfied
         if(Sum > 0){
             res_uncor = CPBayes_uncor( variantName, traitNames, X, s.e., updateDE, MinSlabVar, MaxSlabVar, RP, burn.in )   
             indi_uncor = 1 
         }
     }

     combined_res <- 0
     uncor_use <- 0
     
     if(indi_uncor == 0){ 
          combined_res <- res_cor; uncor_use <- "No"; 
     }else{ 
          combined_res <- res_uncor; uncor_use <- "Yes";
     }

     ptm2 <- proc.time()
     ptm <- ptm2-ptm1
     #print(ptm2-ptm1)
     
     combined_res$uncor_use <- uncor_use                       ## Did the combined strategy choose uncorrelated version
     combined_res$runtime <- ptm
     #combined_res$corCPBayes <- res_cor                        ## collecting the output generated by Correlated CPBayes

     return(combined_res)
 }











