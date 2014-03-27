fRGInteraction <-
function(phenoData,x_A,x_B)
{    
    sample_num=dim(x_A)[1]
    basis_num_A=dim(x_A)[2]
    basis_num_B=dim(x_B)[2]
            
     gamma<-matrix(0,sample_num,basis_num_A*basis_num_B)

     for(i_ in 1:sample_num)
        for(j_ in 1:basis_num_A)
            for(k_ in 1:basis_num_B)
                     gamma[i_,((j_-1)*basis_num_B)+k_]<-x_A[i_,j_]*x_B[i_,k_]

     W<-cbind(x_A,x_B,gamma)
     WTW=t(W)%*%W
     iWTW=ginv(WTW)
     b_hat=iWTW%*%t(W)%*%phenoData
              
     tta=1/(length(phenoData)-basis_num_A-basis_num_B-basis_num_A*basis_num_B)
     delta=tta*as.numeric( t(phenoData)%*%phenoData-t(phenoData)%*%W%*%b_hat )
     varb_hat=delta*iWTW

     gamma_hat<-b_hat[(basis_num_A+basis_num_B+1):length(b_hat),]
     startpos<-basis_num_A+basis_num_B+1
     alpha_hat<-varb_hat[startpos:dim(varb_hat)[1],startpos:dim(varb_hat)[2]]
     Tes=as.numeric(t(gamma_hat)%*%ginv(alpha_hat)%*%gamma_hat)

     rk=qr(alpha_hat)$rank
     rlt=pchisq(Tes,rk,lower.tail=FALSE)

     rlt

}
