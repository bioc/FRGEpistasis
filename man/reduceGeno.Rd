\name{reduceGeno}
\alias{reduceGeno}

\title{
Reduction Dimension of Genotype
}
\description{
Reduce Dimension of Genotype using Functional Regression Model
}
\usage{
reduceGeno(wDir, pheno, gnoFiles, mapFiles, gLst, rng)
}

\arguments{
  \item{wDir}{
The dataset directory. If the dataset is in the working directory, wDir is ".".
}
  \item{pheno}{
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
  \item{rng}{
A numeric value which represents gene region extensible scope.
}
}
\details{

This function reduces the dimension of genotypes of all the genes with Fourier expansion.
In real dataset the genotypes on different chromosome are always organized into different files. And each genotype file is very large.
This function processes the genotype files in turns. The reduced genotype data (the expansion data) of all the chromosomes are combined together.
The expansion data and other information are organized into a list.
During Fourier expansion, the physical position information of the SNPs are used. This is one of merits of our method.
}
\value{
Return a list that includes reduced genotype data, gene names,chromosome information, start index and end index of each gene.
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
gInfo=reduceGeno(wDir,phenoInfo,gnoFiles,mapFiles,gLst,rng)

}
\keyword{Reduction Dimension}
\keyword{Fourier Expansion}
