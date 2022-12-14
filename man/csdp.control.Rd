\name{csdp.control}
\alias{csdp.control}
\title{
Pass control parameters to csdp solver.
}
\description{
Utility function to pass control parameters to csdp solver.
}
\usage{
csdp.control(axtol = 1e-08, 
  atytol = 1e-08, 
  objtol = 1e-08, 
  pinftol = 1e+08, 
  dinftol = 1e+08, 
  maxiter = 100, 
  minstepfrac = 0.9, 
  maxstepfrac = 0.97, 
  minstepp = 1e-08, 
  minstepd = 1e-08, 
  usexzgap = 1, 
  tweakgap = 0, 
  affine = 0, 
  printlevel = 1, 
  perturbobj = 1, 
  fastmode = 0)
}
\arguments{
  \item{axtol}{Tolerance for primal feasibility.
}
  \item{atytol}{Tolerance for dual feasibility.
}
  \item{objtol}{
Tolerance for relative duality gap.
}
  \item{pinftol}{
Tolerance for primal infeasibility.
}
  \item{dinftol}{
Tolerance for dual infeasibility.
}
  \item{maxiter}{
Maximum number of iterations used.
}
  \item{minstepfrac}{
Minimum distance to edge of feasibility region for step.
}
  \item{maxstepfrac}{
Maximum distance to edge of feasibility region for step.
}
  \item{minstepp}{
Failure is declared if primal line search step size is shorter than this parameter.
}
  \item{minstepd}{
Failure is declared if dual line search step size is shorter that this parameter.
}
  \item{usexzgap}{
If 0, then use objective function duality gap.
}
  \item{tweakgap}{
If 1 (and \code{usexzgap}=0) then "fix" negative duality gaps.
}
  \item{affine}{
If 1, only use affine primal-dual steps and do not use barrier function.
} 
  \item{printlevel}{
If 0, no printing, 1 normal printing, higher values result in more debug printing.
}
  \item{perturbobj}{
Amount of objective permutation used.
}
  \item{fastmode}{
If 1, csdp will be faster but also less accurate.
}
}
\details{
Parameters are fully described in CSDP user guide. \url{https://github.com/coin-or/Csdp/}}
\value{
A list with values for all parameters. Any parameters not passed to function are set to default.
}
\references{
\url{https://github.com/coin-or/Csdp/}
}
\author{
Hector Corrada Bravo, CSDP by Brian Borchers
}
\examples{
params <- csdp.control(axtol=1e-6)
}
