### Code contributed by Luciano Selzer

cld <- function (object, ...) {
  UseMethod("cld")
}

cld.lsm <- function(obj, reversed = FALSE, omitNotSig  = FALSE){
  require(plyr)
  require(multcompView)
  require(reshape2)
  #Extract pValues and actual Values
  pValues <- extract_p(obj)
  values <- extract_values.lsm(obj)
  if(is.list(pValues)) {
    # I need to do this to order the levels of according to their mean
    # Otherwise the letters end up mixed
    pValues <- lapply(seq_along(pValues), 
                      function(i) order_p(pValues[[i]], values[[i]]))
    letters <- laply(pValues, wrapLetters, reversed, omitNotSig)
    letters <- melt(letters)
    # Make the order of the levels match those of values data.frame
    letters$Var2 <- factor(letters$Var2, levels = levels(values[[1]][,1]))
    # Order again, otherwise the letters do not match the actual values
    letters <- letters[order(letters[,1], letters[,2]), 3]
  }else {
    pvalues <- order_p(pValues, values)
    letters <- wrapLetters(pValues, reversed = reversed, 
                           omitNotSig = omitNotSig)
  }
  # old behaviour of returning a data.frame
  ans <- cbind(obj[[1]], letters)
  ans
}

extract_p.lsm <- function(x){
  x <- x[[2]]
  ans <- x[,5]
  names(ans) <- gsub(" ", "", row.names(x))
  if(any(grepl("\\|", names(ans)))) {
    s <- do.call(rbind, strsplit(names(ans), split= "\\|"))
    names(ans) <- s[,1]
    # I need the split factor ordered by appearence in the original data.frame
    # It would mix levels of the spliting factor if left to default. 
    # Althougth not everytime.
    f <- factor(s[,2], levels = unique(s[,2]))
    ans <- split(ans, f)
  }
  ans
}

extract_values.lsm <- function(x){
  x <- x[[1]]
  levels(x[, 1]) <- gsub(" ", "", levels(x[,1]))
  # Check if there's a split factor and split if it's there
  if(which(names(x) == "lsmean") == 3) {
    ans <- split(x, x[,2])
  }else ans <- x
  ans
}

order_p <- function(pValues, x){
  # Order the pValues according to their level's mean
  oz <- order(x[,"lsmean"], decreasing = FALSE)
  Lvls <- levels(x[[1]])[oz]
  value <- vec2mat(pValues)
  value <- value[Lvls, Lvls]
  value
}

wrapLetters <- function(x, reversed = FALSE, omitNotSig = FALSE){
  x <- multcompLetters(x, reversed = reversed)[[1]]
  if(omitNotSig){
    if(all(x == x[1])) {
      x[seq_along(x)] <- rep("", length(x))
      
    }
  }
  # Make the order of all the results the same.
  # Because lapply would just rbind the vector together regardless the name
  x <- x[order(names(x))]
  x
}