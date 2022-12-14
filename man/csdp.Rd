\name{csdp}
\alias{csdp}
\alias{csdp_minimal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Solve semidefinite program with CSDP}
\description{
	Interface to CSDP semidefinite programming library. The general statement of the primal problem is
	\deqn{\max\, \mathrm{tr}(CX)}{max tr(CX)}
	\deqn{\mathrm{s.t.}\; A(X) = b}{s.t. A(X) = b}
	\deqn{X \succeq 0}{X >= 0}
	with
	\eqn{A(X)_i = \mathrm{tr}(A_iX)}{A(X)_i = tr(A_iX)}
	where \eqn{X \succeq 0}{X >= 0} means X is positive
	semidefinite, \eqn{C} and all \eqn{A_i} are symmetric matrices of the same
	size and \eqn{b} is a
	vector of length \eqn{m}. 
	
	The dual of the problem is
	\deqn{\min\, b'y}{min b'y}
	\deqn{\mathrm{s.t.}\; A'(y) - C = Z}{s.t. A'(y) - C = Z}
	\deqn{Z \succeq 0}{Z >= 0}

	where
	\eqn{A'(y) = \sum_{i=1}^m y_i A_i.}{A'(y) = \sum_{i=1}^m y_i A_i.}

	Matrices \eqn{C} and \eqn{A_i} are assumed to be block diagonal
	structured, and must be specified that way (see Details). 
}
\usage{
csdp(C, A, b, K,control=csdp.control())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{C}{A list defining the block diagonal cost matrix \eqn{C}.}
  \item{A}{A list of length \eqn{m} containing block diagonal constraint matrices
    \eqn{A_i}. Each constraint matrix \eqn{A_i} is specified by a list
	of blocks as explained in the Details section.}
  \item{b}{A numeric vector of length \eqn{m} containing the right hand side of the constraints.}
  \item{K}{Describes the domain of each block of the sdp problem. It is a list with the following elements:
    \describe{
      \item{type:}{A character vector with entries \code{"s"} or
      \code{"l"} indicating the type of each block. If the \code{j}th entry
      is \code{"s"}, then the \code{j}th block is a positive semidefinite matrix.
      otherwise, it is a vector with non-negative entries.}
      \item{size:}{A vector of integers indicating the dimension of each
      block.}      
      }}
  \item{control}{Control parameters passed to csdp. See CSDP documentation.}
}
\details{
  All problem matrices are assumed to be of block diagonal structure,
  and must be specified as follows:
  \enumerate{
    \item If there are \code{nblocks} blocks specified by \code{K}, then
    the matrix must be a list with \code{nblocks} components.
    \item If \code{K$type == "s"} then the \code{j}th element of the list must define a symmetric
	matrix of size \code{K$size}. It can be an object of class
	\code{"matrix"}, \code{"simple_triplet_sym_matrix"}, or a valid
	class from the class hierarchy in the \code{"Matrix"} package.
   \item If \code{K$type == "l"} then the \code{j}th element of the list
	must be a numeric vector of length \code{K$size}.
  }
  
  This function checks that the blocks in arguments \code{C} and \code{A} agree with
	the sizes given in argument \code{K}. It also checks that the
	lengths of arguments \code{b} and \code{A}
	are equal. It does not check for symmetry in the problem data.
	
	\code{csdp_minimal} is a minimal wrapper to the C code underlying \code{csdp}. It assumes that the arguments 
	\code{sum.block.sizes}, \code{nconstraints}, \code{nblocks}, \code{block.types}, and \code{block.sizes}
	are provided as if they were created by \code{Rcsdp:::prob.info} and that the arguments \code{C}, \code{A}, and 
	\code{b} are provided as if they were created by \code{Rcsdp:::prepare.data}. This function may be useful when 
	calling the csdp functionality iteratively and most of the optimization details stays the same. For example, when the
	control file created by \code{Rcsdp:::write.control.file} stays the same across iterations, but it would be recreated
	on each iteration by \code{csdp}.
}
\value{
  \item{X}{Optimal primal solution \eqn{X}. A list containing blocks in the
	same structure as explained above. Each element is of class
	\code{"matrix"} or a numeric vector as appropriate.}
  \item{Z}{Optimal dual solution \eqn{Z}. A list containing blocks in the same
	structure as explained above. Each element is of class
	\code{"matrix"}
	or a numeric vector as appropriate.}
  \item{y}{Optimal dual solution \eqn{y}. A vector of the same length as
    argument \code{b}}
  \item{pobj}{Optimal primal objective value}
  \item{dobj}{Optimal dual objective value}
  \item{status}{Status of returned solution.
    \describe{
      \item{0:}{Success. Problem solved to full accuracy}
      \item{1:}{Success. Problem is primal infeasible}
      \item{2:}{Success. Problem is dual infeasible}
      \item{3:}{Partial Success. Solution found but full accuracy was
	not achieved}
      \item{4:}{Failure. Maximum number of iterations reached}
      \item{5:}{Failure. Stuck at edge of primal feasibility}
      \item{6:}{Failure. Stuch at edge of dual infeasibility}
      \item{7:}{Failure. Lack of progress}
      \item{8:}{Failure. \eqn{X} or \eqn{Z} (or Newton system \eqn{O}) is singular}
      \item{9:}{Failure. Detected NaN or Inf values}
    }
    }
}
\references{
  \itemize{
  \item \url{https://github.com/coin-or/Csdp/}
  \item Borchers, B.:\cr
      \emph{CSDP, A C Library for Semidefinite Programming} Optimization Methods and Software  11(1):613-623, 1999\cr
      \url{http://euler.nmt.edu/~brian/csdppaper.pdf}
    \item Lu, F., Lin, Y., and Wahba, G.:\cr
      \emph{Robust Manifold Unfolding with Kernel Regularization} TR
      1108, October, 2005. \cr
      \url{http://pages.stat.wisc.edu/~wahba/ftp1/tr1108rr.pdf}
}
}
\author{Hector Corrada Bravo. CSDP written by Brian Borchers. }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
  C <- list(matrix(c(2,1,
                     1,2),2,2,byrow=TRUE),
            matrix(c(3,0,1,
                     0,2,0,
                     1,0,3),3,3,byrow=TRUE),
            c(0,0))
A <- list(list(matrix(c(3,1,
                        1,3),2,2,byrow=TRUE),
               matrix(0,3,3),
               c(1,0)),
          list(matrix(0,2,2),
               matrix(c(3,0,1,
                        0,4,0,
                        1,0,5),3,3,byrow=TRUE),
               c(0,1)))

  b <- c(1,2)
  K <- list(type=c("s","s","l"),size=c(2,3,2))
  csdp(C,A,b,K)

# Manifold Unrolling broken stick example
# using simple triplet symmetric matrices
X <- matrix(c(-1,-1,
              0,0,
              1,-1),nc=2,byrow=TRUE);
d <- as.vector(dist(X)^2);
d <- d[-2]

C <- list(.simple_triplet_diag_sym_matrix(1,3))
A <- list(list(simple_triplet_sym_matrix(i=c(1,2,2),j=c(1,1,2),v=c(1,-1,1),n=3)),
          list(simple_triplet_sym_matrix(i=c(2,3,3),j=c(2,2,3),v=c(1,-1,1),n=3)),
          list(matrix(1,3,3)))

K <- list(type="s",size=3)
csdp(C,A,c(d,0),K)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}

