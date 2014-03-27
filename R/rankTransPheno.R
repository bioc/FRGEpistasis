rankTransPheno <-
function(pheno,para_c)
{
    pheno<-qnorm((rank(pheno)-para_c)/(length(pheno)-2*para_c+1))
    return(pheno)
}
