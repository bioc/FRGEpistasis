\name{outerEpi}
\alias{outerEpi}

\title{
Epistasis Detection Outer Chromosomes
}
\description{
Detect epistasis between 2 genes (or genomic regions) which are on different chromosomes. 
}
\usage{
outerEpi(pheno, gstnd, gStndp, geno_expn, gname, gNamep, gchr, gChrp)
}

\arguments{
  \item{pheno}{
Vector of phenotype data which can be quantitative trait or binary trait.
}
  \item{gstnd}{
Vector of indexes which indicates the start indexes and end indexes of expansed genotype of each gene on one chromosome in the matrix "geno_expn".
}
  \item{gStndp}{
Vector of indexes which indicates the start indexes and end indexes of expansed genotype of each gene on the other chromosome in the matrix "geno_expn".
}
  \item{geno_expn}{
Matrix of expansed genotype data of all the genes.
}
  \item{gname}{
Vector of gene names on one chromosome.
}
  \item{gNamep}{
Vector of gene names on the other chromosome.
}
  \item{gchr}{
Vector of Chromosome number of one chromosome.
}
  \item{gChrp}{
Vector of Chromosome number of the other chromosome.
}
}
\details{
This function tests the epistasis between 2 genes which are on different chromosomes. 
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

smp_num=1000
number_basis<-40
pheno<-sample(c(0:500),smp_num,replace=TRUE)
gname<-c("g1","g2")
gNamep<-c("r1","r2","r3")
gstnd<-c(0,5,9)
gStndp<-c(16,23,29,36)
smpl<-runif(number_basis*smp_num, 0.0, 1.0)
geno_expn<-matrix(smpl,smp_num,number_basis)
gchr<-c(1,1)
gchrp<-c(3,3,3)
outerEpi(pheno,gstnd,gStndp,geno_expn,gname,gNamep,gchr,gchrp)

}

\keyword{ Outer Test  }
\keyword{ Epistasis }
