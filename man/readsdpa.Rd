\name{readsdpa}
\alias{readsdpa}
\alias{writesdpa}
\alias{readsdpa.sol}
\alias{writesdpa.sol}
\title{Reading and writing semidefinite programs for SDPA format files.}
\description{Functions to read and write semidefinite program data and
  solutions in SDPA format.}
\usage{
   readsdpa(file="",verbose=FALSE)
   writesdpa(C,A,b,K,file="")
   readsdpa.sol(K,C,m,file="")
   writesdpa.sol(X,Z,y,K,file="")
}
\arguments{
  \item{file}{The name of the file to read from or write to.}
  \item{C}{Block structured cost matrix}
  \item{A}{List of block structured constraint matrices}
  \item{b}{RHS vector}
  \item{K}{Cone specification, as used in \code{\link{csdp}}}
  \item{X}{Block structured primal optimal solution matrix}
  \item{Z}{Block strucutured dual optimal solution matrix}
  \item{y}{Dual optimal solution vector}
  \item{verbose}{Printout information as problem is read. Passed to CSDP's readsdpa function. Default \code{FALSE}}
  \item{m}{Number of constraints in problem.}
}
\details{
  Block structured matrices must be specified as described in
  \code{\link{csdp}}.
  Files read must be in SDPA format (see
  \url{http://infohost.nmt.edu/~sdplib/FORMAT}).
  However, these functions don't support comments or grouping characters
  (e.g. braces, parentheses) in the block sizes specification.
}
\value{
  Function \code{readsdpa} returns a list with elements \code{C,A,b,K}.
  Function \code{readsdpa.sol} returns a listh with elements
  \code{X,Z,y}.
  All returned matrices are lists of objects of class \code{simple_triplet_sym_matrix}.
}
\references{
  \url{http://infohost.nmt.edu/~sdplib/FORMAT}
}
\author{Hector Corrada Bravo}
\seealso{\code{\link{csdp}}}
\examples{
#  TO DO
}