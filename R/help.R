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

  PPNA <- input$locFDR
  x <- PPNA
  count <- 0
  while(x < 1){ x <- 10*x; count <- count+1 }
  PPNA <- round(PPNA, digits = count+1)

  locFDR <- paste("locFDR", PPNA, sep = " : ")
  cat(locFDR, "\n")
  cat("subset:", " ")
  if(length(input$subset) != 0)
    cat(input$subset, "\n")
  else
    cat("None\n")
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






