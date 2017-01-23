
 ##=================================****** uncorrelated version of CPBayes ******================================##
 ## This function is the main MCMC function for implementing CPBayes in case of uncorrelated summary statistics.
 ## Uncorrelated summary statistics arise for separate case-control statistics without any overlapping subjects. 
 ## It calls the function: initiate_MCMC() from 'common_MCMC_functions.R' to initialize the parameters in the MCMC.
 ## It calls: uncorrelated_beta_update(), Z_update(), q_update(), de_update1(), de_update0() 
 ## from 'uncor_MCMC_functions.R' to update different parameters in MCMC.
 ## It calls select_subset(), overall_pleio_measure() from 'summary_functions.R' to summarize the MCMC data.
 ##=================================**********************************************===============================## 

 #library("MASS")

 CPBayes_uncor = function( variantName, traitNames, X, s.e., updateDE, RP, burn.in )
 {
     ptm1 <- proc.time()
     set.seed(10)

     K = length(X)
     PPAj_thr = 0.25                                           ## PPAj threshold

     tau <- 0.01                                               ## choice of spike sd (var = tau^2)
     nonNullVar <- 1                                           ## central choice of slab variance 
     de <- sqrt(tau^2/nonNullVar)                              ## 1/de = ratio of slab sd and spike sd

     ## v0 (minimum of slab variance - min_var), v1 (maximum of slab variance - max_var)
     min_var <- 0.8
     max_var <- 1.2
     max_de <- tau/sqrt(min_var)                               ## maximum value of 'de'
     min_de <- tau/sqrt(max_var)                               ## mimimum value of 'de'
     shape1_de <- 1                                            ## shape1 parameter of the Beta prior of 'de' (shape2 parameter = 1, always)
     shape1 <- 1                                               ## shape1 parameter for Beta prior of q
     shape2 <- 1                                               ## shape2 parameter for Beta prior of q
     
     ## an informed initialization of the MCMC parameters 
     initiate <- initiate_MCMC( K, X, s.e. )
     beta <- initiate$beta
     Z <- initiate$Z
     q <- initiate$q

     thinning <- 1                                             ## thinning period in the MCMC

     mcmc.samplesize <- (RP-burn.in) %/% thinning; Z.data <- matrix(0,mcmc.samplesize,K); row <- 0 ;
     sim.beta <- matrix(0,mcmc.samplesize,K); sample_probZ_zero <- matrix(0,mcmc.samplesize,K);

     for( rp in 1:RP )
     {
         ## Update beta
         beta = uncorrelated_beta_update(K, s.e., tau, de, Z, X)

         ## Update Z using the q-included version
         res_Z <- Z_update(K, q, tau, de, beta)
         
         ## using the q-integrated out version                                                                                                      
         #res_Z <- Z_integrated_update(K, log_ratio_p1, tau_const, de, beta)     

         ## collecting the otput from Z-update function
         Z <- res_Z$Z
         probZ_zero <- res_Z$prob

         q = q_update(K, Z, shape1, shape2)                    ## Update q
         ## Update de
         if(updateDE == TRUE){
             if(sum(Z) > 0) 
               de <- de_update1(min_de, max_de, shape1_de, beta, Z, tau) 
             else de <- de_update0(min_de, max_de, shape1_de)
         }

         ## collecting the MCMC sample obtained after the burn in period
         if(rp > burn.in && rp%%thinning == 0)
         {   
             row = row + 1 ;
             Z.data[row,] = Z ;
             sim.beta[row,] = beta ;
             sample_probZ_zero[row,] = probZ_zero;
         }
     }                                                        ## closing the loop for MCMC iterations     

     ##----------------------- Compute the summary obtained from the MCMC data ------------------------------------##
     ## selection of subset
     uncor.subset <- select_subset( K, Z.data, mcmc.samplesize )
     selected_traits <- NULL
     if( length(uncor.subset) > 0 ) 
       selected_traits <- traitNames[uncor.subset]
     
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
     log10_BF_uncor <- pleio_evidence$log10_BF
     PPNA.uncor <- pleio_evidence$PPNA    
     
     ptm2 <- proc.time()                                       ## time taken for the analysis
     ptm <- ptm2-ptm1
     #cat(" run time (in seconds):", "\n")
     #print(ptm2-ptm1)

     ## return the outputs. A post summary from the MCMC data can be computed for interesting variants
     data = list( variantName = variantName, log10_BF = log10_BF_uncor, PPNA = PPNA.uncor, 
            subset = selected_traits, important_traits = important_phenos, auxi_data = list( traitNames = traitNames, 
            K = K, mcmc.samplesize = mcmc.samplesize, PPAj = asso.pr, Z.data = Z.data, sim.beta = sim.beta), runtime = ptm )

     return(data)
 }










