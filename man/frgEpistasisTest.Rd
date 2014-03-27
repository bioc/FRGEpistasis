\name{frgEpistasisTest}
\alias{frgEpistasisTest}
\title{
Epistasis Test by Functional Regression Model
}
\description{
 This function is used to analyse the epistasis of genomic region A and genomic region B.
}
\usage{
frgEpistasisTest(pheno, geno_A, pos_A, geno_B, pos_B)
}
\arguments{
  \item{pheno}{
A vector of phenotype which can be quantitative trait or binary trait.
}
  \item{geno_A}{
Genotype matrix of gene ( or genomic region) A.
}
  \item{pos_A}{
Vector of physical positions of SNPs in gene ( or genomic region) A.
}
  \item{geno_B}{
Genotype matrix of gene ( or genomic region) B.
}
  \item{pos_B}{
Vector of physical positions of SNPs in gene ( or genomic region) B.
}
}
\details{
This function is independent with other functions in this package. It is designed for small dataset test.
It takes phenotype, genotype and Physical positions as the input. 
If the position information is NULL, this function considers the SNPs in this gene to be uniformly filled in the gene scope.
First this function expanses the genotypes of gene A and gene B. Then it analyses their epistasis.
}
\value{
 It returns the p value of chi-square test for epistasis detection between gene A and gene B.
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

frgEpistasisTest(pheno,geno_A,pos_A=NULL,geno_B,pos_B=NULL)

}

\keyword{ Epistasis}
\keyword{ Physical Position }
