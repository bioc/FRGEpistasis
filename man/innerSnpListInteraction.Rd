\name{innerSnpListInteraction}
\alias{innerSnpListInteraction}
\title{
Pairwise Interaction Test Inner The Same SNP List
}
\description{
Test the SNP-SNP interaction.
And the SNPs are organized into one data structure.
}
\usage{
innerSnpListInteraction(pheno, snpList)
}

\arguments{
  \item{pheno}{
Vector of phenotype data.
}
  \item{snpList}{
Matrix of the genotypes of all the SNPs for testing the pairwise interactions.
}
}
\details{
This function aims to test the pairwise interactions between the SNPs  organized into the same data structure.
It takes phenotype and genotypes of the SNPs as the input. 
And output all the p values for the interactions of SNP pairs.

}
\value{
Return a frame contains names of all the SNPs pairs and  p values for interactions of these pairs.
}

\author{
Futao Zhang
}

\examples{

pheno<- round(runif(1000,40,60))
geno<- as.data.frame(matrix(round(runif(5000,0,2)),1000,5))
innerSnpListInteraction(pheno,geno)

}

\keyword{ SNP-SNP Interaction}
\keyword{ Pairwise Interaction }
