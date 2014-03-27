\name{fRGInteraction}
\alias{fRGInteraction}
\title{
Interaction Test by Functional Regression Model
}
\description{
Test interaction between two gene (or genomic regions) with chi-squared test.
}
\usage{
fRGInteraction(phenoData, x_A, x_B)
}

\arguments{
  \item{phenoData}{
Vector of phenotype data which can be quantitative trait or binary trait.
}
  \item{x_A}{
Expansion data matrix of Genotype of gene A. 
}
  \item{x_B}{
Expansion data matrix of Genotype of gene B.
}
}
\details{
This function takes phenotype vector and expansed genotype matrices as input.
It is the most important part of this software package. It is called by functions "innerEpi" and "innerEpi" of this package.
The interaction  between gene A and gene B is tested with chi-squared test.
}
\value{
It returns the p value of chi-squared test for epistasis detection between gene A and gene B.
}

\author{
Futao Zhang
}

\examples{
x_A<-as.matrix(rnorm(1000,mean=0,sd=1))
x_B<-as.matrix(rnorm(1000,mean=0,sd=1))
phenoData<-runif(1000,15,60)
fRGInteraction(phenoData,x_A,x_B)
}

\keyword{central chi-square distribution }

