\name{fourierExpansion}
\alias{fourierExpansion}

\title{
Fourier Expansion of Genotype
}
\description{
This function aims to expand the genotype of one gene (or genomic region)
with Fourier Expansion.
}
\usage{
fourierExpansion(gene_idx, geno, gene_list, snp_map, rng)
}

\arguments{
  \item{gene_idx}{
The expansion gene index in the gene annotation list.
}
  \item{geno}{
Genotype of all the genes in the gene annotation list.
}
  \item{gene_list}{
Gene annotation list which includes gene name, chromosome, start position and end position.
}
  \item{snp_map}{
SNP genetic map includes chromosome, snp indentifier, genetic distance and base-pair position.
}
  \item{rng}{
A numeric value which represents gene region extensible scope.
}
}
\details{
This function reduces the dimension of one gene(or genomic region) with Fourier Expansion. 
Fist extract out the genotype of this gene with the gene annotation and the SNP map information.
Then expanse the gene with the genotypes and SNP positions if the number of SNPs in the gene is over 3.
Otherwise the raw genotypes of the gene would be returned.
The number of Fourier Basis is selected to explain 80 percent of genetic variation.
}
\value{
If the SNPs number of the gene is over 3, returns the expansion of the genotype which is a matrix 
with the dimension Sample number * Fourier Basis number.
If the SNPs number of the gene is no more than 3, returns the raw genotype of the gene which is a matrix 
with the dimension Sample number * SNP number.
}

\author{
Futao Zhang
}

\examples{
gLst<-read.csv(system.file("extdata", "gene.list.csv", package="FRGEpistasis"))
fdata<-read.table(system.file("extdata", "simGeno-chr1.raw", package="FRGEpistasis"),header=TRUE)
geno<-fdata[,-1:-6]
snp_map<-read.table(system.file("extdata", "chr1.map", package="FRGEpistasis"))
fourierExpansion(1, geno, gLst, snp_map, 0)
}

\keyword{ Fourier Expansion }
\keyword{ Basis Number }
