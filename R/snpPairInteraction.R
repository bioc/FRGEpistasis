snpPairInteraction <-
function(pheno,snp1,snp2)
{
    spl_num=length(pheno)
    vec=rep(0,spl_num)
    for(k in 1:spl_num)
    {
        vec[k]= snp1[k]*snp2[k];
    }
    W<-cbind(rep(1,spl_num),snp1,snp2,vec)
    WTW=t(W)%*%W
    iWTW=ginv(WTW)
    beta_hat=iWTW%*%t(W)%*%pheno
    Y_betahat <-pheno-W%*%beta_hat
    delta=as.numeric(t(Y_betahat )%*%Y_betahat)/(spl_num-4)
    var_beta_hat=delta*iWTW
    T_sta=beta_hat[4]/sqrt(var_beta_hat[4,4])
    pval=pnorm(T_sta, mean=0, sd=1, lower.tail=FALSE, log.p=FALSE)
    return(pval)
}
