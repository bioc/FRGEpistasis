pointwiseInteraction <-
function(phenoData,x_A,x_B)
{    
    idx=is.na(phenoData)
        phenoData=phenoData[!idx]
    x_A=as.matrix(x_A[!idx,])
     x_B=as.matrix(x_B[!idx,])
         
     spl_num=dim(x_A)[1]           
        k1=dim(x_A)[2]
    k2=dim(x_B)[2]
    vec=rep(0,spl_num)
    min_pval=1

        for(i_ in 1:k1)
    {
        for(j_ in 1:k2)
        {
            for(k_ in 1:spl_num)
            {
                vec[k_]= x_A[k_,i_]*x_B[k_,j_];
            }
            W<-cbind(rep(1,spl_num),x_A[,i_],x_B[,j_],vec)
            WTW=t(W)%*%W
                 iWTW=ginv(WTW)
                 beta_hat=iWTW%*%t(W)%*%phenoData
            Y_betahat <-phenoData-W%*%beta_hat
            delta=as.numeric(t(Y_betahat )%*%Y_betahat)/(spl_num-4)
                 var_beta_hat=delta*iWTW
            tst_sta=beta_hat[4]/sqrt(var_beta_hat[4,4])
            if(!is.na(tst_sta))
            {
                pval=pnorm(tst_sta,mean=0,sd=1,lower.tail=FALSE,log.p=FALSE)
                if(pval<min_pval) min_pval=pval
            }            
        }
    }    
    return(min_pval)
}
