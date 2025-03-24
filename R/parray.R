##' Parse a conditional formula
##'
##' Parse a formula of the form `y ~ d_1+d_2+...+d_n | c_1+c_2+...+c_m`
##' to extract the response variable `y`, dependent variables `d_i`,
##' and conditioning variables `c_i` as character vectors.
##'
##' @title Parse a conditional formula
##' @param formula A formula
##' @return A list with elements rvar, dvars, and cvars
parse_cformula <- function(formula) {
  # Ensure formula is of class "formula"
  if (!inherits(formula, "formula")) stop("Input must be a formula")

  # Extract the left-hand side (lhs) and right-hand side (rhs)
  lhs <- if(length(formula) > 2) formula[[2]]
  rhs <- formula[[length(formula)]]

  # Extract the weight (lhs of ~)
  rvar <- if(!is.null(lhs)) all.vars(lhs)

  # Check if `|` is present in the rhs
  if(inherits(rhs, "call") && identical(as.character(rhs[[1]]),"|")) {
    dvars <- all.vars(rhs[[2]])
    cvars <- all.vars(rhs[[3]])
    if(any(cvars %in% dvars))
      warning("Invalid formula: dependent and conditioning variables must be disjoint")
  } else {
    dvars <- all.vars(rhs)
    cvars <- NULL
  }

  # Return as a named list
  list(rvar = rvar, dvars = dvars, cvars = cvars)
}


##' Create a named array of probabilities from a data frame or named
##' array using a formula interface.
##'
##' Given a formula of the form `w ~ d_1+d_2+...+d_n | c_1+c_2+...+c_m`,
##' `ptabs` creates a named array of probabilities of the dependent
##' variables `d_i` conditioned on the conditioning variables `c_i`.
##' If data is a dataframe, the optional response `w` specifies a
##' weighting variable, otherwise it is ignored.  If data is a
##' parray object, the result is conditioned on both the `c_i` and
##' any existing conditioning variables in `data`.
##'
##' Where conditioning results in division by zero, the undefined
##' values are replaced with the elements of `fill`, which may be a
##' vector or another named array.
##'
##' @title Create a named array of probabilities
##' @param formula A formula
##' @param data A data frame or named array
##' @param fill A vector or named array of replacement values for
##'   division by zero.
##' @importFrom stats reformulate xtabs
##' @return A named array of probabilities
##' @export
ptabs <- function(formula,data,fill=0) {
  ## Extract variables
  vlist <- parse_cformula(formula)
  ## Create table if data is a data frame
  if(inherits(data,"data.frame")) {
    data <- xtabs(reformulate(c(vlist$dvars,vlist$cvars),response=vlist$rvar),data)
  } else {
    if(!is.null(vlist$rvar)) warning("response variable is ignored for named arrays")
  }
  ## Check data is a named array
  if(!is.array(data) && is.null(dimnames(data)))
    stop("data must be a named array or data frame")
  ## Check variables from formula are in data
  vars <- names(dimnames(data))
  fvars <- c(vlist$dvars,vlist$cvars)
  if(length(miss <- setdiff(fvars,vars)) > 0L)
    stop("Variables not found in data: ",paste(miss,collapse=" "))
  ## Check any existing conditioning
  cdims <- max(0,attr(data,"cdims"))
  ddims <- length(vars)-cdims
  cvars <- vars[ddims+seq_len(cdims)]
  if(any(vlist$dvars %in% cvars)) stop("Data is conditioned on dependent variables")
  ## Marginalize if necessary
  if(length(setdiff(vars[seq_along(ddims)],fvars)) > 0L)
    data <- marginalize(data,fvars)
  ## Convert to array of conditional probabilities
  condition(data,vlist$cvars,fill=fill)
}

##' Condition an array of weights or probabilities on a given
##' set of variables.
##'
##' This function conditions an array of weights or probabilities on a
##' given set of variables.  The array is permuted so that the
##' conditioning variables are the first dimensions of the array, and
##' then attempts to normalize so that for each combination of the
##' conditioning variables the corresponding array slice sums to one.
##' Where this normalization results in division by zero, the
##' corresponding values are replaced with the element of "fill",
##' which may be a vector or named array.
##'
##' @title Conditioning probability arrays
##' @param P A named array of weights or probabilities
##' @param by A vector of variable names to condition on
##' @param fill A vector or named array of replacement values for
##'   division by zero.
##' @return An array of conditional probabilities
##' @export
condition <- function(P,by=c(),fill=0) {
  ## Get variable names
  vars <- names(dimnames(P))
  ## Indices of conditioning variables
  idx <- match(by,vars)
  cdims <- max(0,attr(P,"cdims"))
  ddims <- length(vars)-cdims
  if(cdims > 0L) idx <- union(ddims+seq_len(cdims),idx)
  cdims <- length(idx)
  ## Reorder and normalize
  P <- aperm(P,c(setdiff(seq_along(vars),idx),idx))
  if(cdims==0) {
    S <- sum(P)
  } else {
    dm <- dim(P)
    ddims <- length(dm)-cdims
    S <- colSums(P,dims=ddims)
    S <- matrix(S,nrow=prod(dm[seq_len(ddims)]),ncol=length(S),byrow=TRUE)
  }
  P <- P/as.vector(S)
  ## Fill in div by zero
  if(any(zero <- (S==0))) {
   if(inherits(fill,"parray")) {
    ## Check fill array is compatible
    vars <- names(dimnames(P))
    fvars <- names(dimnames(fill))
    if(length(setdiff(fvars,vars)) > 0L)
      stop("fill array has incompatible variables")
    ## Extend fill array
    ext <- match(setdiff(vars,fvars),vars)
    fill <- aperm(array(fill,c(dim(fill),dim(P)[ext])),match(vars,c(fvars,vars[ext])))
    ## Fill adjusting for replication of dependent fill variables
    scale <- prod(dim(P)[ext[ext<=length(vars)-cdims]])
    P[zero] <- fill[zero]/scale
   } else {
    P[zero] <- fill
   }
  }
  ## Record conditioning variables
  attr(P,"cdims") <- length(idx)
  class(P) <- c("parray","array")
  P
}


##' Marginalize an array of weights or probabilities to a given
##' set of variables.
##'
##' Marginalize an array of weights or probabilities, summing over
##' any dependent variables not given in `to`.  Any conditioning
##' variables are retained.  If `prob` is `TRUE`, the result is
##' normalized so that for each combination of the conditioning
##' variables, the sum over the variables in `to` is one.
##'
##' @title Marginalizing probability arrays
##' @param P A named array of weights or probabilities
##' @param to A vector of variable names to marginalize to
##' @return An array of marginal probabilities
##' @export
marginalize <- function(P,to) {
  ## Indices that are retained
  dnms <- dimnames(P)
  idx <- match(to,names(dnms))
  cdims <- max(0,attr(P,"cdims"))
  ddims <- length(dnms)-cdims
  if(any(idx > ddims)) warning("Marginalizing conditioning variables")
  idx <- c(idx,ddims+seq_len(cdims))
  P <- as.array(apply(P,idx,sum))
  dimnames(P) <- dnms[idx]
  attr(P,"cdims") <- cdims
  class(P) <- c("parray","array")
  P
}


##' Compute the product of two arrays of conditional probabilities or
##' weights.
##'
##' Given P1 = P(X|Y,Z) and P2 = P(Y|Z), this function computes the
##' product P(X,Y|Z) = P(X|Y,Z) P(Y|Z).
##'
##' @title Product of two arrays
##' @param P1 A parray of weights or probabilities
##' @param P2 A parray of weights or probabilities
##' @return A parray.
##' @export
product <- function(P1,P2) {
  ## Get variable names
  dnms1 <- dimnames(P1)
  dnms2 <- dimnames(P2)
  vars1 <- names(dnms1)
  dims1 <- dim(P1)
  vars2 <- names(dnms2)
  dims2 <- dim(P2)

  ## Check dimensions
  vars <- intersect(vars1,vars2)
  if(any(dims1[match(vars,vars1)]!=dims2[match(vars,vars2)]))
    stop("Incompatible dimensions")

  ## Dependent and conditioning variables
  cdims1 <- max(0,attr(P1,"cdims"))
  ddims1 <- length(vars1)-cdims1
  cvars1 <- vars1[ddims1+seq_len(cdims1)]
  dvars1 <- vars1[seq_len(ddims1)]
  cdims2 <- max(0,attr(P2,"cdims"))
  ddims2 <- length(vars2)-cdims2
  cvars2 <- vars2[ddims2+seq_len(cdims2)]
  dvars2 <- vars2[seq_len(ddims2)]

  if(length(intersect(dvars1,dvars2)) > 0L)
    stop("Incompatible dependent variables")

  ## Variables in the product
  dvars <- c(dvars1,dvars2)
  cvars <- c(cvars2,setdiff(cvars1,dvars2))
  vars <- c(dvars,cvars)
  dims <- c(dims1,dims2)[match(vars,c(vars1,vars2))]

  ## Extend P1 and P2
  ext <- match(setdiff(vars,vars1),vars)
  P1 <- aperm(array(P1,dim=c(dims1,dims[ext])),match(vars,c(vars1,vars[ext])))
  ext <- match(setdiff(vars,vars2),vars)
  P2 <- aperm(array(P2,dim=c(dims2,dims[ext])),match(vars,c(vars2,vars[ext])))
  P <- P1*P2
  dimnames(P) <- c(dnms1,dnms2)[vars]
  attr(P,"cdims") <- length(cvars)
  class(P) <- c("parray","array")
  P
}


##' Select a subarray from a named array, indexing by variable names.
##'
##' @title Named indexing
##' @param narray A named array
##' @param ... Variable values to extract
##' @return A subarray
##' @export
subarray <- function(narray, ...) {
  args <- list(...)
  sub <- dimnames(narray)
  sub[names(args)] <- args
  do.call(`[`,c(list(narray),sub))
}


##' Convert a parray to a data frame.
##'
##' Convert a parray to a data frame.
##'
##' @title Convert a parray to a data frame
##' @param x A parray
##' @param row.names The names of the rows
##' @param optional ignored
##' @param responseName The name of the variable
##' @param ... currently ignored
##' @return A data frame
##' @export
as.data.frame.parray <- function(x,row.names=NULL,optional,responseName="(Weight)",...) {
  df <- do.call(expand.grid,lapply(dimnames(x),function(v) factor(v,levels=v)))
  df[[responseName]] <- as.vector(x)
  rownames(df) <- row.names
  df
}


