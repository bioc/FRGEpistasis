\name{rankTransPheno}
\alias{rankTransPheno}

\title{
Rank-Based Inverse Normal Transformation
}
\description{
Rank-Based Inverse Normal Transformation of Phenotype
}
\usage{
rankTransPheno(pheno, para_c)
}

\arguments{
  \item{pheno}{
Vector of phenotype which is the quantitative trait.
}
  \item{para_c}{
Adjust parameter, commonly as 0,1/3,3/8 or 1/2.
}
}
\details{
Some variables are not normally distributed.
And using statistical tests on this data can give misleading results because they do not meet the statistical assumptions. 
This function implements Rank-Based Inverse Normal Transformation to make phenotype normally distributed.
}
\value{
Return vector of rank-based inverse normal transformed phenotype.
}
\references{
T. Mark Beasley, Stephen Erickson and David B. Allison. Rank-Based Inverse Normal Transformations are Increasingly Used, But are They Merited?
Behav Genet. 2009 Sep.;39(5):580-595.
}
\author{
Futao Zhang
}

\examples{

c=0.5
smp_num=100
pheno<-sample(c(0:500),smp_num,replace=TRUE)
rankTransPheno(pheno,c)

}

\keyword{ Rank-Based}
\keyword{ Inverse Normal Transformation}
