\name{innerEpi}
\alias{innerEpi}
\title{
Epistasis Detection Inner one Chromosome
}
\description{
Detect epistasis between 2 genes (or genomic regions) both of which are on the same chromosome. 
}
\usage{
innerEpi(pheno, gstnd, geno_expn, gname, gchr)
}

\arguments{
  \item{pheno}{
Vector of phenotype data which can be quantitative trait or binary trait.
}
  \item{gstnd}{
Vector of indexes which indicates the start indexes and end indexes of expansed genotype of each gene on current chromosome in the matrix "geno_expn".

}
  \item{geno_expn}{
Matrix of expansed genotype data of all the genes.
}
  \item{gname}{
Vector of gene names on current chromosome.
}
  \item{gchr}{
Vector of Chromosome number of current chromosome.
}
}
\details{
This function tests the epistasis between 2 genes both of which are on the same chromosome. 
It takes expansed genotype data as the input. 
First the data of the gene are extracted from "geno_expn" with "gstnd" and "gname".
Then the function "fRGInteraction" will be called.
}
\value{
Return a matrix which contains the gene names of the gene pairs and the p values of chi-squared test for the epistasis of the gene pairs.
}

\author{
Futao Zhang
}

\examples{

smp_num<-1000
number_basis<-9
pheno<-sample(c(0:500),smp_num,replace=TRUE)
gname<-c("g1","g2")
gstnd<-c(0,5,9)
smpl<-runif(number_basis*smp_num, 0.0, 1.0)
geno_expn<-matrix(smpl,smp_num,number_basis)
gchr<-c(1,1)
innerEpi(pheno,gstnd,geno_expn,gname,gchr)

}

\keyword{ Inner Test }
\keyword{ Epistasis}
