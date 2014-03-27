\name{snpPairInteraction}
\alias{snpPairInteraction}

\title{
SNP-SNP interaction
}
\description{
Test the interaction of one SNP with another
}
\usage{
snpPairInteraction(pheno, snp1, snp2)
}

\arguments{
  \item{pheno}{
Vector of phenotype data which can be quantitative trait or binary trait.
}
  \item{snp1}{
Vector of genotype data of SNP1.
}
  \item{snp2}{
Vector of genotype data of SNP2.
}
}
\details{
This function tests the interaction of one SNP with another. 
}
\value{
Return the p value for snp-snp interaction
}

\author{
Futao Zhang
}

\examples{

pheno<- round(runif(1000,40,60))
snp1<-round(runif(1000,0,2))
snp2<-round(runif(1000,0,2))
pval=snpPairInteraction(pheno,snp1,snp2)

}

\keyword{ SNP-SNP interaction}

