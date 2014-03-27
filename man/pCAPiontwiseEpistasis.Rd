\name{pCAPiontwiseEpistasis}
\alias{pCAPiontwiseEpistasis}

\title{
Epistasis Test by PCA Method and Piontwise Method
}
\description{
This function is the another entrance of the software package.
It tests the Genome Wide Epistasis by PCA Method and Piontwise Method.
}
\usage{
pCAPiontwiseEpistasis(wDir, oEpi, phenoInfo, gnoFiles, mapFiles, gLst, rng)
}

\arguments{
  \item{wDir}{
The dataset directory. If the dataset is in the working directory, wDir is ".".
}
  \item{oEpi}{
Output data frame which contains all the names of gene pairs and the p values for their epistasis.
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
  \item{rng}{
A numeric value which represents gene region extensible scope.
}
}
\details{

The genotypes on different chromosome are stored in different files. 
The full names of these files are listed in the index file that is taken as the input parameter.
After the index file is loaded, the function knows where to read the genotypes files.
This function analyses the genotypes files one by another. 
That means this function tests the epistasis of genomic regions both of which are on the same chromosome(file),
then the epistasis of genomic regions which are on different chromosomes(files).
This function can test epistasis both with PCA method and pointwise method. 
For a pair of genes, we assume that the total number of all possible SNP pairs is K, 
The minum p value for SNP-SNP interaction among the K pairs is output as the pointwise method result of the gene pair. 
}

\value{
Return a data frame which contains all the names of gene pairs and the p values of chi-square test for their epistasis.
}

\author{
Futao Zhang
}

\examples{

work_dir <-paste(system.file("extdata", package="FRGEpistasis"),"/",sep="")
##read the list of genotype files
geno_files<-read.table(system.file("extdata", "list_geno.txt", package="FRGEpistasis"))
##read the list of map files
mapFiles<-read.table(system.file("extdata", "list_map.txt", package="FRGEpistasis"))
##read the phenotype file
phenoInfo <- read.csv(system.file("extdata", "phenotype.csv", package="FRGEpistasis"),header=TRUE)
##read the gene annotation file
gLst<-read.csv(system.file("extdata", "gene.list.csv", package="FRGEpistasis"))
##define the extension scope of gene region
rng=0

##log transformation
phenoInfo [,2]=log(phenoInfo[,2])
out_epi<-data.frame()

pCAPiontwiseEpistasis(work_dir,out_epi,phenoInfo,geno_files,mapFiles,gLst,rng)

}
\keyword{ Epistasis}
\keyword{ PCA}
\keyword{ Pointwise interaction}

