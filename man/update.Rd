%    Copyright (c) 2012-2016 Russell V. Lenth                                %

\name{update}
\alias{update}
\alias{update.ref.grid}
\alias{lsm.options}
\alias{get.lsm.option}
\alias{defaults.lsm}
\alias{pmm.options}
\alias{get.pmm.option}

\title{
Set or retrieve options for objects and summaries in \pkg{lsmeans}
}
\description{
Objects of class \code{ref.grid} or \code{lsmobj} contain several settings in their \code{"misc"} slot that affect primarily the 
defaults used by \code{\link{summary}}. This \code{update} method allows them to be changed more safely than by modifying this slot directly. 

In addition, the user may set defaults for all objects using \samp{options(lsmeans = ...)}, or more conveniently using the \code{lsm.options} and \code{get.lsm.option} functions documented here (or its courtesy wrappers, \code{pmm.options} and \code{get.pmm.option} for those who dislike the \sQuote{least-squares means} terminology).
}
\usage{
\S3method{update}{ref.grid}(object, ..., silent = FALSE)

lsm.options(...)
get.lsm.option(x, default = defaults.lsm[[x]])

pmm.options(...)
get.pmm.option(...)
}
\arguments{
  \item{object}{An object of class \code{ref.grid} (or its extension, \code{lsmobj})
}
  \item{\dots}{
Arguments specifying elements' names and their new values.
}
  \item{silent}{If \code{FALSE}, a message is displayed for any unmatched names.}
  \item{x}{Character string holding an option name for \code{lsm.options}.}
  \item{default}{Return value if \code{x} is not found.}
}

\details{

\subsection{Using \code{update}}{
In \code{update}, the names in \code{\dots} are partially matched against those that are valid, and if a match is found, it adds or replaces the current setting. The valid names are

\describe{
\item{\code{tran}, \code{tran2}}{(\code{list} or \code{character}) specifies the transformation which, when inverted, determines the results displayed by \code{\link{summary}}, \code{\link{predict}}, or \code{\link{lsmip}} when \code{type="response"}. The value may be the name of a standard transformation from \code{\link{make.link}} or additional ones supported by name, such as \code{log2}; or, for a custom transformation, a \code{list} containing at least the functions \code{linkinv} (the inverse of the transformation) and \code{mu.eta} (the derivative thereof). The \code{\link{make.tran}} function returns such lists for a number of popular transformations. See the help page of \code{\link{make.tran}} for details as well as information on the additional named transformations that are supported. \code{tran2} is just like \code{tran} except it is a second transformation (i.e., a response transformation in a generalized linear model).}

\item{\code{tran.mult}}{Multiple for \code{tran}. For example, for the response transformation \samp{2*sqrt(y)} (or \samp{sqrt(y) + sqrt(y + 1)}, for that matter), we should have \code{tran = "sqrt"} and \code{tran.mult = 2}. If absent, a multiple of 1 is assumed.}

\item{\code{estName}}{(\code{character}) is the column label used for displaying predictions or LS means.}

\item{\code{inv.lbl}}{(\code{character)}) is the column label to use for predictions or LS means when \code{type="response"}.}

\item{\code{by.vars}}{(\code{character)} vector or \code{NULL}) the variables used for grouping in the summary, and also for defining subfamilies in a call to \code{\link{contrast}}.}

\item{\code{pri.vars}}{(\code{character} vector) are the names of the grid variables that are not in \code{by.vars}. Thus, the combinations of their levels are used as columns in each table produced by \code{\link{summary}}.}

\item{\code{alpha}}{(numeric) is the default significance level for tests, in \code{\link{summary}} as well as \code{\link{cld}} and \code{\link{plot}} when \samp{intervals = TRUE}}

\item{\code{adjust}}{(\code{character)}) is the default for the \code{adjust} argument in \code{\link{summary}}.}

\item{\code{estType}}{(\code{character}) is the type of the estimate. It should match one of \samp{c("prediction","contrast","pairs")}. This is used along with \code{"adjust"} to determine appropriate adjustments to P values and confidence intervals.}

\item{\code{famSize}}{(integer) is the \code{nmeans} parameter for \code{\link{ptukey}} when \code{adjust="tukey"}. }

\item{\code{infer}}{(\code{logical} vector of length 2) is the default value of \code{infer} in \code{\link{summary}}.}

\item{\code{level}}{(numeric) is the default confidence level, \code{level}, in \code{\link{summary}}}

\item{\code{df}}{(numeric) overrides the default degrees of freedom with a specified single value.}

\item{\code{null}}{(numeric) null hypothesis for \code{summary} or \code{test} (taken to be zero if missing).}

\item{\code{side}}{(numeric or character) \code{side} specification for for \code{summary} or \code{test} (taken to be zero if missing).}

\item{\code{delta}}{(numeric) \code{delta} specification for \code{summary} or \code{test} (taken to be zero if missing).}

\item{\code{predict.type}}{(character) sets the default method of displaying predictions in \code{\link{summary}}, \code{\link{predict}}, and \code{\link{lsmip}}. Valid values are \code{"link"} (with synonyms \code{"lp"} and \code{"linear"}), or \code{"response"}.}

\item{\code{avgd.over}}{(\code{character)} vector) are the names of the variables whose levels are averaged over in obtaining marginal averages of predictions, i.e., LS means. Changing this might produce a misleading printout, but setting it to \code{character(0)} will suppress the \dQuote{averaged over} message in the summary.}

\item{\code{initMesg}}{(\code{character}) is a string that is added to the beginning of any annotations that appear below the \code{\link{summary}} display.}

\item{\code{methDesc}}{(\code{character}) is a string that may be used for creating names for a list of \code{lsmobj} objects. }

\item{(any slot name)}{If the name matches an element of \code{slotNames(object)}, that slot is replaced by the supplied value, if it is of the required class (otherwise an error occurs). Note that all the other possibilities above refer to elements of \code{misc}; hence, you probably don't want to replace \code{misc} itself. The user must be very careful in replacing slots because they are interrelated; for example, the \code{levels} and \code{grid} slots must involve the same variable names, and the lengths and dimensions of \code{grid}, \code{linfct}, \code{bhat}, and \code{V} must conform. 
}
} % end \describe
} % end subsection

\subsection{Using \code{lsm.options}}{
In \code{lsm.options}, we may set or change the default values for the above attributes in the \code{lsmeans} option list(see \code{\link{options}}). Currently, the following elements of this list are used if specified:
\describe{
\item{\code{ref.grid}}{A named \code{list} of defaults for objects created by \code{\link{ref.grid}}. This could affect other objects as well. For example, if \code{lsmeans} is called with a fitted model object, it calls \code{ref.grid} and this option will affect the resulting \code{lsmobj} object.}
\item{\code{lsmeans}}{A named \code{list} of defaults for objects created by \code{\link{lsmeans}} (or \code{\link{lstrends}}).}
\item{\code{contrast}}{A named \code{list} of defaults for objects created by \code{\link{contrast}} (or \code{\link{pairs}}).}
\item{\code{summary}}{A named \code{list} of defaults used by the methods \code{\link{summary}}, \code{\link{predict}}, and \code{\link{lsmip}}. The only option that can affect the latter two is \code{"predict.method"}.}
\item{\code{estble.tol}}{Tolerance for determining estimability in rank-deficient cases. If absent, the value in \code{defaults.lsm$estble.tol)} is used.}
\item{(others)}{Other options may be accessed by support code for particular model classes (see \code{\link{models}}). For example, the \code{lmer.df}, \code{disable.pbkrtest}, and \code{pbkrtest.limit} options affect how degrees of freedom are computed for \code{lmerMod} objects (\pkg{lme4} package).}
} % end \describe
} % end subsection

} % end details
\value{
\code{update} returns a copy of \code{object} with its \code{"misc"} slot modified (and perhaps other slots). \code{lsm.options} returns the current options (same as the result of \samp{getOption("lsmeans")}) -- invisibly, unless called with no arguments.
}


\author{
Russell V. Lenth
}
\note{If a call to \code{\link{lsmeans}}, \code{\link{contrast}}, or \code{\link{ref.grid}} contains a non-\code{NULL} \code{options} list, those options are passed in a call to \code{update} on the constructed object before it is returned. This allows you, for example, to override the defaults used by \code{\link{summary}}. In addition, user defaults may be set using an \code{link{options}} setting for \code{"lsmeans"}. It should be a list with one or more named elements \code{lsmeans}, \code{contrast}, or \code{ref.grid}, used for setting the defaults for objects constructed by functions of these same names. Note that options can get \dQuote{inherited}. See the examples.

Unlike the \code{update} method for model classes (\code{lm}, \code{glm}, etc.), this does not re-fit or re-estimate anything; but it does affect how \code{object} is treated by other methods for its class.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{summary}}, \code{\link{make.tran}}
}
\examples{
# An altered log transformation

warp.lm1 <- lm(log(breaks + 1) ~ wool*tension, data = warpbreaks)
rg1 <- update(ref.grid(warp.lm1), 
        tran = list(linkinv = function(eta) exp(eta) - 1,
                    mu.eta = function(eta) exp(eta)),
        inv.lbl = "pred.breaks")
        
summary(rg1, type = "response")

\dontrun{
lsm.options(ref.grid = list(level = .90),
            contrast = list(infer = c(TRUE,FALSE)),
            estble.tol = 1e-6)
# Sets default confidence level to .90 for objects created by ref.grid
# AS WELL AS lsmeans called with a model object (since it creates a 
# reference grid). In addition, when we call 'contrast', 'pairs', etc.,
# confidence intervals rather than tests are displayed by default.
}

\dontrun{
lsm.options(disable.pbkrtest = TRUE)
# This forces use of asymptotic methods for lmerMod objects.
# Set to FALSE or NULL to re-enable using pbkrtest.
}

# See tolerance being used for determining estimability
get.lsm.option("estble.tol")
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }% __ONLY ONE__ keyword per line
