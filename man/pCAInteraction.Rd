\name{pCAInteraction}
\alias{pCAInteraction}

\title{
Epistasis Test by Principal Component Analysis 
}
\description{
Test the epistasis between two genes (or genomic regions) with the principal components analysis method.
}
\usage{
pCAInteraction(phenoData, x_A, x_B)
}

\arguments{
  \item{phenoData}{
Vector of phenotype data which can be quantitative trait or binary trait.
}
  \item{x_A}{
 Matrix of genotype of gene A.
}
  \item{x_B}{
 Matrix of genotype of gene B.
}
}
\details{
This function takes phenotype vector and genotype matrices as input and tests the epistasis using PCA method. 
The number of principal components is determined by PCA to explain 80 percent of the genetic variation.
The interaction  between gene A and gene B is tested with chi-squared test.
}
\value{
It returns the p value of chi-squared test for epistasis detection between gene A and gene B.
}

\author{
Futao Zhang
}

\examples{

smp_num=1000
number_snp_A=25
number_snp_B=20
pheno<-sample(c(0:500),smp_num,replace=TRUE)
smpl=rep(0,number_snp_A*smp_num)
idx_1=sample(c(1:(number_snp_A*smp_num)),ceiling(number_snp_A*smp_num/100))
idx_2=sample(c(1:(number_snp_A*smp_num)),ceiling(number_snp_A*smp_num/200))
smpl[idx_1]=1
smpl[idx_2]=2
geno_A=matrix(smpl,smp_num,number_snp_A)

smpl=sample(c(0,1,2),number_snp_B*smp_num,replace=TRUE)
geno_B=matrix(smpl,smp_num,number_snp_B)
pCAInteraction(pheno,geno_A,geno_B)
}

\keyword{ Epistasis }
\keyword{ PCA }
