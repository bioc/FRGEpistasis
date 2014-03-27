\name{outerSnpListInteraction}
\alias{outerSnpListInteraction}

\title{
Pairwise Interaction Test Outer The SNP Lists
}
\description{
Test the SNP-SNP interaction.
And the SNPs are organized into two different SNP Lists.
}
\usage{
outerSnpListInteraction(pheno, snpList1, snpList2)
}

\arguments{
  \item{pheno}{
Vector of phenotype data.
}
  \item{snpList1}{
Matrix of the genotypes of all the SNPs on the first SNP list for testing the pairwise interactions.
}
  \item{snpList2}{
Matrix of the genotypes of all the SNPs on the second SNP list for testing the pairwise interactions.
}
}
\details{
This function aims to test the pairwise interactions between the SNPs  organized into different data structures.
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
snp_list_1 <- as.data.frame(matrix(round(runif(3000,0,2)),1000,3))
snp_list_2 <- as.data.frame(matrix(round(runif(5000,0,2)),1000,5))
colnames(snp_list_1 )<-c("rs10","rs11","rs12")
colnames(snp_list_2 )<-c("rs20","rs21","rs22","rs23","rs24")
pheno<- round(runif(1000,40,60))
outerSnpListInteraction(pheno,snp_list_1,snp_list_2)
}

\keyword{ SNP-SNP interaction  }
\keyword{ Pairwise Interaction }
