##============================ Some supporting functions for functions for users ===========================##

## Check primary variables
checkPrimaryVar <-  function(VAR, nameVAR)
{
  # Check whether VAR is a vector
  if(!is.vector(VAR))
  {
    if(is.matrix(VAR) && any(dim(VAR)==1))
    {
      VAR <- as.vector(VAR)
      warning( paste(nameVAR, "is a matrix!"), call. = FALSE)
    }
    else
      stop(paste(nameVAR, "must be a vector."), call. = FALSE)
  }
  # Check whether VAR is a numeric vector
  if(!is.numeric(VAR))
    stop(paste(nameVAR, "must be a numeric vector."), call. = FALSE)
  # Check whether there is any NA
  if(any(is.na(VAR)))
    stop(paste(nameVAR, "for one or more phenotypes are missing!"), call. = FALSE)
  # Check whether there is more than one non-missing arguments
  if(length(VAR) <= 1)
    stop(paste("Number of elements in the", nameVAR, "vector must be more than 1!"), call. = FALSE)
  return(VAR)
}

## Check input variable 'Phenotypes'
checkPhen <-  function(Phenotypes, BetaHat)
{
  # Check whether argument 3 is a vector
  if(!is.vector(Phenotypes))
    stop("Phenotypes must be a vector.", call. = FALSE)
  # Check whether argument 3 is a character vector
  if(!is.character(Phenotypes))
    stop("Phenotypes must be a character vector.", call. = FALSE)
  # Check whether there is duplicate phenotyes
  if(length(Phenotypes) > length(unique(Phenotypes)))
    stop("Two or more phenotypes have the same name!", call. = FALSE)
  # Check whether argument 3 is a vector of length more than 1
  if(length(Phenotypes) != length(BetaHat))
    stop("BetaHat and Phenotypes vectors must have the same number of elements!", call. = FALSE)
}

## Check input variable
checkVarName <- function(variantName)
{
  # Check whether argument 4 is a vector
  if(!is.vector(variantName))
    stop("Variant must be a vector.", call. = FALSE)
  # Check whether argument 4 is a vector of length 1
  if(length(variantName) > 1)
    stop("Variant must be a vector of length 1.", call. = FALSE)
  # Check whether argument 4 is NA
  if(is.na(variantName))
  {
    warning("Variant is NA!", call. = FALSE)
    variantName <- as.character(variantName)
  }
  # Check whether argument 4 is not NA but numeric
  if(!is.na(variantName) && !is.character(variantName))
  {
    warning("Variant is not a character vactor!", call. = FALSE)
    variantName <- as.character(variantName)
  }
    return(variantName)
}

checkCorln <- function(Corln, BetaHat)
{
  # Check whether Corln is a data.frame
  if(is.data.frame(Corln))
    stop("Corln must be a matrix not a data.frame. Use as.matrix() to convert the data.frame into matrix.", call. = FALSE)
  # Check whether Corln is a matrix
  if(!is.matrix(Corln))
    stop("Corln must be a matrix.", call. = FALSE)
  # Check whether Corln is numeric
  if(!is.numeric(Corln))
    stop("Corln must be a numeric matrix.", call. = FALSE)
  # Check whether there is any
  if(any(is.na(Corln)))
    stop("One or more entries of Corln are missing!", call. = FALSE)
  # Check for Corln, whether number of rows = number of columns
  if(nrow(Corln) != ncol(Corln))
    stop("Number of rows and columns of Corln are different!", call. = FALSE)
  # Check whether number of rows of corln matrix is same as no. of entries in BetaHat
  if(nrow(Corln) != length(BetaHat))
    stop("Number of rows of Corln and length of BetaHat do not match!", call. = FALSE)
  # Save as matrix
    Cor <- as.matrix(Corln)
  # Check whether a symmetric matrix. First, make row names and col names identical.
    row.names(Cor) <- colnames(Cor)
  if(!isSymmetric(Cor))
    stop("Corln is not symmetric!", call. = FALSE)
  # Check whether a negative definite matrix
  if(det(Cor) < 0)
    stop("Corln is negative definite!", call. = FALSE)
  # Check whether diagonal elements are 1
  if(dist(rbind(diag(Cor), rep(1, dim(Cor)[1]))) != 0)
    stop("Diagonal elements of Corln are not 1!", call. = FALSE)
  # Check whether a singular matrix
  if(det(Cor) == 0)
    warning("Corln is a singular matrix!", call. = FALSE)

  return(Cor)
}

print_result <- function(input)
{
  #cat("RESULT ::", "\n")
  #gvar <- paste(" genetic_variant", input$genetic_variant, sep = " : ")
  #cat(gvar, "\t")
  #BF <- input$log10_BF
  #BF <- round(BF, digits = 2)
  #BF <- paste("log10_BF", BF, sep = ": ")
  #cat(BF, "\n")

  #PPNA <- input$locFDR
  #x <- PPNA
  #count <- 0
  #while(x < 1){ x <- 10*x; count <- count+1 }
  #PPNA <- round(PPNA, digits = count+1)

  #locFDR <- paste("locFDR", PPNA, sep = " : ")
  #cat(locFDR, "\n")
  cat("Important traits with trait-specific posterior prob of assoc: \n")
  if(sum(dim(input$important_traits)) > 0){
    dat = input$important_traits
    dat$PPAj = round(dat$PPAj, 2)
    print(dat)
    cat("\n")
  }else{
    cat("None\n")
  }
}


chk_n_stp1 <-  function(arg, argName)
{
  # Check whether arg is a data.frame
  if(is.data.frame(arg))
    stop(paste(argName, "must be a matrix not a data.frame. Use as.matrix() to convert the data.frame into matrix."), call. = FALSE)
  # Check whether arg is a matrix
  if(!is.matrix(arg))
    stop(paste(argName, "must be a matrix."), call. = FALSE)
  # Check for arg, whether number of rows = number of columns
  if(nrow(arg) != ncol(arg))
    stop(paste("Number of rows and columns of", argName, "are different!"), call. = FALSE)
  # Check whether there is any
  if(any(is.na(arg)))
    stop(paste("One or more entries of", argName, "are missing!"), call. = FALSE)
  # Check whether arg is numeric
  if(!is.numeric(arg))
    stop(paste(argName, "must be a numeric matrix."), call. = FALSE)
  # Check whether non-negative integer matrix
  if(any(arg%%1 != 0) || any(arg < 0))
    stop(paste("Every element of", argName, "must be a non-negative integer."), call. = FALSE)

}

chkEstCorln <- function(n11, n00, n10)
{
  chk_n_stp1(n11, "n11")
  chk_n_stp1(n00, "n00")
  chk_n_stp1(n10, "n10")
  if((nrow(n11) != nrow(n00)) || (nrow(n10) != nrow(n00)))
    stop("n11, n00, n10 matrices must have same dimension.", call. = FALSE)
  row.names(n11) <- colnames(n11)
  row.names(n00) <- colnames(n00)
  row.names(n10) <- colnames(n10)
  if(!isSymmetric(n11))
    stop("n11 must be symmetric!", call. = FALSE)
  if(!isSymmetric(n00))
    stop("n00 must be symmetric!", call. = FALSE)
  if(any(diag(n11) == 0))
    stop("Diagonal elements of n11 must be positive integer.", call. = FALSE)
  if(any(diag(n00) == 0))
    stop("Diagonal elements of n00 must be positive integer.", call. = FALSE)
  if(any(diag(n10) != 0))
    stop("Diagonal elements of n10 must be zero.", call. = FALSE)
}

chkSlabVarSpikeVar <- function(SpikeVar=0.0001, SlabVar=0.8){
  defaultVal <- list("SpikeVar" = 0.0001, "SlabVar" = 0.8)
  minVal <- list("SpikeVar" = 0, "SlabVar" = 0)
  var_list <- list("SpikeVar" = defaultVal[["SpikeVar"]], "SlabVar" = defaultVal[["SlabVar"]])
  if(!missing(SpikeVar)) var_list[["SpikeVar"]] <- SpikeVar
  if(!missing(SlabVar)) var_list[["SlabVar"]] <- SlabVar
  for(var in names(var_list)) {
    if(!is.vector(var_list[[var]]) || length(var_list[[var]]) != 1) {
      warning(var, " is not a scalar (default option used).", call. = FALSE)
      var_list[[var]] <- defaultVal[[var]]
    }
    else if(!is.numeric(var_list[[var]])) {
      warning(var, " is not numeric (default option used).", call. = FALSE)
      var_list[[var]] <- defaultVal[[var]]
    }
    else if(var_list[[var]] <= minVal[[var]]) {
      warning(var, " is not positive (default option used).", call. = FALSE)
      var_list[[var]] <- defaultVal[[var]]
    }
  }
  return(var_list)
}

chkSlabVar <- function(MinSlabVar = 0.6, MaxSlabVar = 1.0) {
  defaultVal <- list("MinSlabVar" = 0.6, "MaxSlabVar" = 1.0)
  boundVal <- list("MinSlabVar" = 0.01, "MaxSlabVar" = 10.0)
  var_list <- list("MinSlabVar" = defaultVal[["MinSlabVar"]], "MaxSlabVar" = defaultVal[["MaxSlabVar"]])
  if(!missing(MinSlabVar))  var_list[["MinSlabVar"]] <- MinSlabVar
  if(!missing(MaxSlabVar))  var_list[["MaxSlabVar"]] <- MaxSlabVar
  for(var in names(var_list)) {
    # Check whether argument 6 is a vector of length 1
    if(!is.vector(var_list[[var]]) || (length(var_list[[var]]) != 1)) {
      warning(var, " is not a scalar (default option used).", call. = FALSE)
      var_list[[var]] <- defaultVal[[var]]
    }
    # Check whether argument 6 is numeric
    else if(!is.numeric(var_list[[var]])) {
      warning(var, " is not numeric (default option used).", call. = FALSE)
      var_list[[var]] <- defaultVal[[var]]
    }
  }
  # Check whether argument 6 is more than its minimum bound
  if(var_list[["MinSlabVar"]] < boundVal[["MinSlabVar"]]){
    warning("MinSlabVar should not be smaller than ", boundVal[["MinSlabVar"]], ", so it is assigned to ", boundVal[["MinSlabVar"]], ".", call. = FALSE)
    var_list[["MinSlabVar"]] <- boundVal[["MinSlabVar"]]
  }
  if(var_list[["MaxSlabVar"]] > boundVal[["MaxSlabVar"]]){
    warning("MaxSlabVar should not be bigger than ", boundVal[["MaxSlabVar"]], ", so it is assigned to ", boundVal[["MaxSlabVar"]], ".", call. = FALSE)
    var_list[["MaxSlabVar"]] <- boundVal[["MaxSlabVar"]]
  }
  # Argument 6 and 7 :: Checking MinSlabVar < MaxSlabVar
  if(var_list[["MinSlabVar"]] >= var_list[["MaxSlabVar"]]){
    print(var_list[["MinSlabVar"]])
    print(var_list[["MaxSlabVar"]])
    warning("MaxSlabVar is not bigger than MinSlabVar! (default option used).", call. = FALSE)
    var_list[["MaxSlabVar"]] <- defaultVal[["MaxSlabVar"]]
    var_list[["MinSlabVar"]] <- defaultVal[["MinSlabVar"]]
  }
  return(var_list)
}


chkMCMCparam <- function(MCMCiter = 20000, Burnin = 5000) {
  defaultVal <- list("MCMCiter" = 20000, "Burnin" = 5000)
  lBound <- list("MCMCiter" = 7000, "Burnin" = 2000)
  var_list <- list("MCMCiter" = defaultVal[["MCMCiter"]], "Burnin" = defaultVal[["Burnin"]])
  if(!missing(MCMCiter)) var_list[["MCMCiter"]] <- MCMCiter
  if(!missing(Burnin)) var_list[["Burnin"]] <- Burnin
  for(var in names(var_list)) {
    if(!is.vector(var_list[[var]]) || (length(var_list[[var]]) != 1)){
      warning(var, " is not a scalar (default option used).", call. = FALSE)
      var_list[[var]] <- defaultVal[[var]]
    }else if(!is.numeric(var_list[[var]]) || var_list[[var]]%%1 != 0){
      print(var_list[[var]])
      warning(var, " not provided as integer (default option used).", call. = FALSE)
      var_list[[var]] <- defaultVal[[var]]
    }else if(var_list[[var]] < lBound[[var]]){
      warning(var, " should be at least ", lBound[[var]], " (default option used).", call. = FALSE)
      var_list[[var]] <- defaultVal[[var]]
    }
  }
  if((var_list[["MCMCiter"]] - var_list[["Burnin"]]) < 5000)
  {
    warning("MCMC sample size (MCMCiter - Burnin) provided less than 5000 (default options used)", call. = FALSE)
    var_list[["MCMCiter"]] <- defaultVal[["MCMCiter"]]
    var_list[["Burnin"]] <- defaultVal[["Burnin"]]
  }
  return(var_list)
}






