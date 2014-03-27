\name{fRGEpistasis}
\alias{fRGEpistasis}

\title{
 Genome Wide Epistasis Study by Functional Regress Model
}
\description{
This function is the entrance of the software package.
It tests the Genome Wide Epistasis by Functional Regression Model.
}
\usage{
fRGEpistasis(wDir, phenoInfo, gnoFiles, mapFiles, gLst, fdr, rng)
}

\arguments{
  \item{wDir}{
The dataset directory. If the dataset is in the working directory, wDir is ".".
}
  \item{phenoInfo}{
  It is a matrix with two columns. One column is the individual ID and the other is the phenotype.
  The phenotype can be quantitative trait or binary trait.
}
  \item{gnoFiles}{
The vector of genotype file names. It contains the genotype file names indicating where to read the genotype files.
}
  \item{mapFiles}{
The vector of SNP genetic map file names. It contains the map file names indicating where to read the genetic map files.
}
  \item{gLst}{
Gene annotation which includes gene name, chromosome, start position and end position.
}
  \item{fdr}{
FDR control threshold, When this value == 1, turn FDR control off. 
}
  \item{rng}{
A numeric value which represents gene region extensible scope.
}
}
\details{
Firstly this package reduces the dimension of genotype of all the genomic regions.
Secondly this package tests the epistasis of genomic regions both of which are on the same chromosome(file).
Thirdly this package tests the epistasis of genomic regions which are on different chromosomes(files).

This function is memory efficient with high performance.
Memory efficiency: Only store reduced expansion data of genotypes instead of raw data of genotypes.
    This package reduces the dimension of genotype of all the genomic regions(see details of function "reduceGeno").
	In real dataset the genotypes on different chromosome are always organized into different files. And each genotype file is very large.
	Reading all the files into memory is unacceptable. This package reads the files one by one and reduces the genotype dimension with Fourier expansion.
	In order to inform the package how many files and where to read, we need two data structures "gnoFiles" and "mapFiles" to store the file names.
	
high performance:
	Each data file only needs to read  once and reduce dimension once. So I/O times are reduced and repeated computing of data reduction was avoided.
	This method is a kind of group test. We take a gene(or genomic region) as the test unit. The number of Test is sharply reduced comparing with point-wise interaction (SNP-SNP) test.
    The dimension of genotype is reduced by functional expansion, So the time of each test is reduced.

}
\value{
Return a data frame which contains all the names of gene pairs and the p values of chi-square test for their epistasis.
}

\author{
Futao Zhang
}

\examples{
wDir <-paste(system.file("extdata", package="FRGEpistasis"),"/",sep="")
gnoFiles<-read.table(system.file("extdata", "list_geno.txt", package="FRGEpistasis"))
mapFiles<-read.table(system.file("extdata", "list_map.txt", package="FRGEpistasis"))
phenoInfo <- read.csv(system.file("extdata", "phenotype.csv", package="FRGEpistasis"),header=TRUE)
gLst<-read.csv(system.file("extdata", "gene.list.csv", package="FRGEpistasis"))
rng=0
fdr=0.05
out_epi <- data.frame( )
phenoInfo [,2]=log(phenoInfo [,2])
out_epi = fRGEpistasis(wDir,phenoInfo,gnoFiles,mapFiles,gLst,fdr,rng)
}

\keyword{ Epistasis }
\keyword{ Functional Regress Model}
\keyword{ Dimension Reduction}
