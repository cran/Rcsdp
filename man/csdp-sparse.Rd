\name{csdp-sparse}
\alias{simple_triplet_sym_matrix-class}
\alias{simple_triplet_sym_matrix}
\alias{as.simple_triplet_sym_matrix.matrix}
\alias{as.matrix.simple_triplet_sym_matrix}
\alias{as.vector.simple_triplet_sym_matrix}
\alias{.simple_triplet_zero_sym_matrix}
\alias{.simple_triplet_diag_sym_matrix}
\alias{.simple_triplet_random_sym_matrix}
\title{Simple support for sparse matrices}
\description{Support for sparse matrices in package \code{Rcsdp}. The class
  \code{simple_triplet_sym_matrix} is defined to provide support
  for symmetric sparse matrices. It's definition is copied from the package \code{relations} by Kurt
  Hornik. Coercion functions from objects of class \code{matrix} and
  classes in the \code{Matrix} hierarchy are provided.}
\usage{
   simple_triplet_sym_matrix(i,j,v,n=max(c(i,j)),check.ind=FALSE)
   \method{as.simple_triplet_sym_matrix}{matrix}(x,check.sym=FALSE)
   \method{as.matrix}{simple_triplet_sym_matrix}(x,...)
   \method{as.vector}{simple_triplet_sym_matrix}(x,...)
   .simple_triplet_zero_sym_matrix(n,mode="double")
   .simple_triplet_diag_sym_matrix(x,n)
   .simple_triplet_random_sym_matrix(n,occ=.1,nnz=occ*n*(n+1)/2,rfun=rnorm,seed=NULL,...)
}
\arguments{
  \item{i}{Row indices of non-zero entries.}
  \item{j}{Column indices of non-zero entries.}
  \item{v}{Non-zero entries.}
  \item{n}{Size of matrix.}
  \item{check.ind}{Checks that arguments \code{i} and \code{j} indicate entries in the lower triangular part of the matrix. Default \code{FALSE}.}
  \item{check.sym}{Checks if matrix object is symmetric. Default \code{FALSE}.}
  \item{x}{Object of class \code{matrix} or \code{simple_triplet_sym_matrix}.}
  \item{mode}{Type of zero matrix to create. Default \code{double}.}
  \item{occ}{Ratio of occupancy of random sparse matrix. Default \code{.1}.}
  \item{nnz}{Number of non-zero entries in random sparse matrix. Default corresponds to \code{occ=.1}.}
  \item{rfun}{Function to generate random entries in sparse matrix. Default \code{rnorm}.}
  \item{seed}{Random number generator seed. Set by function \code{set.seed} before generating random sparse matrix. Default \code{NULL}.}
  \item{...}{Arguments passed on to casting functions.}
}
\details{
  TO DO
}
\value{
  TO DO
}
\references{
  TO DO
}
\author{Hector Corrada Bravo}
\seealso{\link{csdp}}
\keyword{utilities}
\examples{
#  TO DO
}