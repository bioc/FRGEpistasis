\name{logTransPheno}
\alias{logTransPheno}

\title{
 logarithmic transformation
}
\description{
Logarithmic Transformation of Phenotype
}
\usage{
logTransPheno(pheno)
}

\arguments{
  \item{pheno}{
Vector of phenotype which is the quantitative trait.
}
}
\details{
Some variables are not normally distributed.
And using statistical tests on this data can give misleading results because they do not meet the statistical assumptions. 
Many variables have log-normal distributions. 
}
\value{
Return vector of transformed phenotype.
}

\examples{

smp_num=100
pheno<-sample(c(0:500),smp_num,replace=TRUE)
logTransPheno(pheno)
}

\keyword{ logarithmic transformations  }

