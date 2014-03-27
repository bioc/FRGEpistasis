\name{pointwiseInteraction}
\alias{pointwiseInteraction}

\title{
Pointwise Interaction Test
}
\description{
Test the epistasis of the gene pair by pointwise method
}
\usage{
pointwiseInteraction(phenoData, x_A, x_B)
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

This function takes phenotype vector and genotype matrices as input and tests the epistasis using pointwise method.
For a pair of genes, we assume that the total number of all possible SNP pairs is K (one SNP from one gene and the other SNP from the other gene).
The interaction of each SNP pair between the two genes is tested.
The minum p value for SNP-SNP interaction among the K pairs is output as the pointwise method result of the gene pair.  

}
\value{
Return the minum p value for SNP-SNP interaction among the K pairs
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
pointwiseInteraction(pheno,geno_A,geno_B)

}
\keyword{ Epistasis }
\keyword{ pointwise interaction}
\keyword{ SNP pair}
