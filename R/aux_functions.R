#Auxiliary functions

metric2info <- function(metric) {
  if (startsWith(metric, "L1")) "L1 statistic"
  else if (startsWith(metric, "L2")) "L2 statistic"
  else if (metric == "mahal") "average pairwise Mahalanobis distance"
  else if (metric == "euclid") "average pairwise Euclidean distance"
  else if (metric == "custom") "average pairwise distance"
  else if (metric == "energy") "energy distance"
}

#Function to turn a vector into a string with "," and "and" or "or" for clean messages. 'and.or'
#controls whether words are separated by "and" or "or"; 'is.are' controls whether the list is
#followed by "is" or "are" (to avoid manually figuring out if plural); quotes controls whether
#quotes should be placed around words in string. From MatchIt.
word_list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately
  L <- length(word.list)
  word.list <- add_quotes(word.list, quotes)

  if (L == 0) {
    out <- ""
    attr(out, "plural") <- FALSE
  }
  else {
    word.list <- word.list[!word.list %in% c(NA_character_, "")]
    L <- length(word.list)
    if (L == 0) {
      out <- ""
      attr(out, "plural") <- FALSE
    }
    else if (L == 1) {
      out <- word.list
      if (is.are) out <- paste(out, "is")
      attr(out, "plural") <- FALSE
    }
    else {
      and.or <- match_arg(and.or)
      if (L == 2) {
        out <- paste(word.list, collapse = paste0(" ", and.or," "))
      }
      else {
        out <- paste(paste(word.list[seq_len(L-1)], collapse = ", "),
                     word.list[L], sep = paste0(", ", and.or," "))

      }
      if (is.are) out <- paste(out, "are")
      attr(out, "plural") <- TRUE
    }

  }
  return(out)
}

#Add quotation marks around a string. From MatchIt.
add_quotes <- function(x, quotes = 2) {
  if (!isFALSE(quotes)) {
    if (isTRUE(quotes) || as.integer(quotes) == 2) x <- paste0("\"", x, "\"")
    else if (as.integer(quotes) == 1) x <- paste0("\'", x, "\'")
    else stop("'quotes' must be boolean, 1, or 2.")
  }
  x
}

#More informative and cleaner version of base::match.arg. From MatchIt.
match_arg <- function(arg, choices, several.ok = FALSE) {
  #Replaces match_arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg))
    stop("No argument was supplied to match_arg.", call. = FALSE)
  arg.name <- paste(deparse(substitute(arg), width.cutoff = 500L), collapse = " ")

  if (missing(choices)) {
    formal.args <- formals(sys.function(sysP <- sys.parent()))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }

  if (is.null(arg))
    return(choices[1L])
  else if (!is.character(arg))
    stop(paste0("The argument to '", arg.name, "' must be NULL or a character vector"), call. = FALSE)
  if (!several.ok) {
    if (identical(arg, choices))
      return(arg[1L])
    if (length(arg) > 1L)
      stop(paste0("The argument to '", arg.name, "' must be of length 1"), call. = FALSE)
  }
  else if (length(arg) == 0)
    stop(paste0("The argument to '", arg.name, "' must be of length >= 1"), call. = FALSE)

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L))
    stop(paste0("The argument to '", arg.name, "' should be ", if (length(choices) > 1) {if (several.ok) "at least one of " else "one of "} else "",
                word_list(choices, and.or = "or", quotes = 2), "."),
         call. = FALSE)
  i <- i[i > 0L]
  if (!several.ok && length(i) > 1)
    stop("There is more than one match in 'match_arg'")
  choices[i]
}

#Turn a vector into a 0/1 vector. 'zero' and 'one' can be supplied to make it clear which is
#which; otherwise, a guess is used. From MatchIt.
binarize <- function(variable, zero = NULL, one = NULL) {
  if (length(unique(variable)) > 2) stop(paste0("Cannot binarize ", paste(deparse(substitute(variable)), collapse = " "), ": more than two levels."))
  if (is.character(variable) || is.factor(variable)) {
    variable <- factor(variable, nmax = 2)
    unique.vals <- levels(variable)
  }
  else {
    unique.vals <- unique(variable, nmax = 2)
  }

  if (is.null(zero)) {
    if (is.null(one)) {
      if (can_str2num(unique.vals)) {
        variable.numeric <- str2num(variable)
      }
      else {
        variable.numeric <- as.numeric(variable)
      }

      if (0 %in% variable.numeric) zero <- 0
      else zero <- min(variable.numeric, na.rm = TRUE)

      return(setNames(as.integer(variable.numeric != zero), names(variable)))
    }
    else {
      if (one %in% unique.vals) return(setNames(as.integer(variable == one), names(variable)))
      else stop("The argument to 'one' is not the name of a level of variable.", call. = FALSE)
    }
  }
  else {
    if (zero %in% unique.vals) return(setNames(as.integer(variable != zero), names(variable)))
    else stop("The argument to 'zero' is not the name of a level of variable.", call. = FALSE)
  }
}

#Determine whether a character vector can be coerced to numeric
can_str2num <- function(x) {
  nas <- is.na(x)
  suppressWarnings(x_num <- as.numeric(as.character(x[!nas])))
  return(!anyNA(x_num))
}

#Cleanly coerces a character vector to numeric; best to use after can_str2num()
str2num <- function(x) {
  nas <- is.na(x)
  suppressWarnings(x_num <- as.numeric(as.character(x)))
  x_num[nas] <- NA
  return(x_num)
}

#Make interaction vector out of matrix of covs. From MatchIt.
exactify <- function(X, nam = NULL, sep = "|", include_vars = FALSE) {
  if (is.null(nam)) nam <- rownames(X)
  if (is.matrix(X)) X <- setNames(lapply(seq_len(ncol(X)), function(i) X[,i]), colnames(X))
  if (!is.list(X)) stop("X must be a matrix, data frame, or list.")

  #Ensure no ambiguity is created by sep
  sep0 <- sep
  unique.x <- unlist(lapply(X, function(x) as.character(unique(x))))
  while (any(grepl(sep, unique.x, fixed = TRUE))) {
    sep0 <- paste0(sep0, sep)
  }

  if (include_vars) {
    for (i in seq_along(X)) {
      if (is.character(X[[i]]) || is.factor(X[[i]])) {
        X[[i]] <- paste0(names(X)[i], ' = "', X[[i]], '"')
      }
      else {
        X[[i]] <- paste0(names(X)[i], ' = ', X[[i]])
      }
    }
  }

  out <- do.call("paste", c(X, sep = sep0))
  if (!is.null(nam)) names(out) <- nam
  out
}

#Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#Get covariates (RHS) vars from formula. From MatchIt.
get.covs.matrix <- function(formula = NULL, data = NULL) {

  if (is.null(formula)) {
    fnames <- colnames(data)
    fnames[!startsWith(fnames, "`")] <- paste0("`", fnames[!startsWith(fnames, "`")], "`")
    formula <- reformulate(fnames)
  }
  else formula <- update(terms(formula, data = data), NULL ~ . + 1)

  mf <- model.frame(terms(formula, data = data), data,
                    na.action = na.pass)

  chars.in.mf <- vapply(mf, is.character, logical(1L))
  mf[chars.in.mf] <- lapply(mf[chars.in.mf], factor)

  X <- model.matrix(formula, data = mf,
                    contrasts.arg = lapply(Filter(is.factor, mf),
                                           contrasts, contrasts = FALSE))
  assign <- attr(X, "assign")[-1]
  X <- X[,-1,drop=FALSE]
  attr(X, "assign") <- assign

  return(X)
}

#Generalized inverse; port of MASS::ginv()
generalized_inverse <- function(sigma) {
  sigmasvd <- svd(sigma)
  pos <- sigmasvd$d > max(1e-8 * sigmasvd$d[1L], 0)
  sigma_inv <- sigmasvd$v[, pos, drop = FALSE] %*% (sigmasvd$d[pos]^-1 * t(sigmasvd$u[, pos, drop = FALSE]))
  return(sigma_inv)
}

#Effective sample size
ESS <- function(w) {
  sum(w)^2/sum(w^2)
}

#Weighted mean, faster than weighted.mean()
w_m <- function(x, w = NULL) {
  if (is.null(w)) sum(x)/length(x)
  else sum(x*w)/sum(w)
}

#Weighted KS statistic, from cobalt's col_w_ks()
w_ks <- function(x, treat, w = NULL) {
  if (is.null(w)) w <- rep(1, length(treat))

  tval1 <- treat[1]

  w_ <- w
  w_[treat == tval1] <-  w[treat == tval1]/sum(w[treat == tval1])
  w_[treat != tval1] <- -w[treat != tval1]/sum(w[treat != tval1])

  ord <- order(x)
  cumv <- abs(cumsum(w_[ord]))[c(diff(x[ord]) != 0, TRUE)]
  ks <- if (length(cumv) == 0) 0 else max(cumv)

  return(ks)
}

#Get calling function
get_calling_function <- function() {
  package.funs <- getNamespaceExports(packageName()) #Note: doesn't capture S3 methods
  function.stack <- unlist(lapply(sys.calls(), function(x) deparse1(x[[1]])))

  matching.functions <- package.funs[package.funs %in% function.stack]

  if (length(matching.functions) == 0) return(character(0))
  else return(paste0(matching.functions[1], "()"))
}

#Wrapper for poly() that doesn't error when degree >= length(unique(x))
poly_ <- function(x, ..., degree = 1, coefs = NULL, raw = FALSE, simple = FALSE) {
  nu <- length(unique(x))
  if (nu <= 2) return(x)
  out <- stats::poly(x, ..., degree = min(degree, nu - 1), coefs = coefs, raw = raw, simple = simple)
  # class(out) <- "matrix"
  out
}

#Process data to be safe for fitting model and getting predictable results without error
process_safe_for_model <- function(data, formula = NULL) {

  if (is.null(formula)) vars <- names(data)
  else vars <- all.vars(formula)

  for (i in vars) {
    nu <- length(unique(data[[i]]))
    if (nu == 1L) data[[i]] <- 1
    else if (nu == 2L) data[[i]] <- as.numeric(data[[i]] == data[[i]][1])
    else if (is.character(data[[i]]) || is.factor(data[[i]])) data[[i]] <- factor(data[[i]])
  }
  data
}

#Safely compute robust SE CIs when there are non-finite entries in meat
#by using simpler variance
safe_ci <- function(fit, treatment, vcov. = sandwich::vcovHC, type = "HC3", alpha = .05, ...) {

  v <- vcov.(fit, type = type, ...)

  while (type != "const") {

    if (is.finite(v[treatment, treatment])) {
      return(lmtest::coefci(fit, treatment, vcov. = v, level = 1 - alpha))
    }

    type <- switch(type,
                   "HC5" = "HC4m",
                   "HC4m" = "HC4",
                   "HC4" = "HC3",
                   "HC3" = "HC2",
                   "HC2" = "HC1",
                   "HC1" = "HC0",
                   "HC0" = "const")

    v <- vcov.(fit, type = type, ...)
  }

  lmtest::coefci(fit, treatment, vcov. = v, level = 1 - alpha)
}