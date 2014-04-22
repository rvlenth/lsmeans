\name{oranges}
\alias{oranges}
\docType{data}
\title{
Orange sales
}
\description{
This example dataset on sales of oranges has two factors, two covariates, and two responses. There is one observation per factor combination.
}
\usage{data(oranges)}
\format{
  A data frame with 36 observations on the following 6 variables.
  \describe{
    \item{\code{store}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6}. The store that was observed.}
    \item{\code{day}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6}. The day the observation was taken (same for each store).}
    \item{\code{price1}}{a numeric vector. Price of variety 1.}
    \item{\code{price2}}{a numeric vector. Price of variety 2.}
    \item{\code{sales1}}{a numeric vector. Sales (per customer) of variety 1.}
    \item{\code{sales2}}{a numeric vector. Sales (per customer) of variety 2.}
  }
}
%\details{}
\source{
Download from \url{http://ftp.sas.com/samples/A56655}.
}
\references{
Littell, R., Stroup W., Freund, R. (2002) \emph{SAS For Linear Models} (4th edition). SAS Institute. ISBN 1-59047-023-0. 
}
\examples{
require(lsmeans)

# Example on p.244 of Littell et al.
oranges.lm <- lm(sales1 ~ price1*day, data = oranges)
lsmeans(oranges.lm, "day")

# Example on p.246
lsmeans(oranges.lm, "day", at = list(price1 = 0))
}
\keyword{datasets}
