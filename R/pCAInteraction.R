pCAInteraction <-
function(phenoData,x_A,x_B)
{    
    cMean_A<-colMeans(x_A)
    for(i in 1:dim(x_A)[2])
        x_A[,i]<-x_A[,i]-cMean_A[i]

    eigen_A<-eigen(var(x_A))
    scores_A<-x_A%*%eigen_A$vectors
     percentage_A=(cumsum(eigen_A$values))/sum(eigen_A$values)    
    a=percentage_A<=0.8
    a[1]=TRUE
     zeta=as.matrix(scores_A[,a])
    
    cMean_B<-colMeans(x_B)
    for(i in 1:dim(x_B)[2])
        x_B[,i]<-x_B[,i]-cMean_B[i]


    eigen_B<-eigen(var(x_B))
    scores_B<-x_B%*%eigen_B$vectors
     percentage_B=(cumsum(eigen_B$values))/sum(eigen_B$values)    
    b=percentage_B<=0.8
    b[1]=TRUE
     eta=as.matrix(scores_B[,b])

    gamma<-matrix(0,dim(zeta)[1],dim(zeta)[2]*dim(eta)[2])

     for(i_ in 1:dim(zeta)[1])
        for(j_ in 1:dim(zeta)[2])
             for(k_ in 1:dim(eta)[2])
                 gamma[i_,((j_-1)*dim(eta)[2])+k_]<-zeta[i_,j_]*eta[i_,k_]

     W<-cbind(zeta,eta,gamma)
     WTW=t(W)%*%W
     iWTW=ginv(WTW)
     b_hat=iWTW%*%t(W)%*%phenoData
     smpSize<-length(phenoData)
     tta<-1/(smpSize-dim(zeta)[2]-dim(eta)[2]-dim(zeta)[2]*dim(eta)[2])
     delta=tta*as.numeric(t(phenoData)%*%phenoData-t(phenoData)%*%W%*%b_hat)
     varb_hat=delta*iWTW

     gamma_hat<-b_hat[(dim(zeta)[2]+dim(eta)[2]+1):length(b_hat),]
     startpos<-dim(zeta)[2]+dim(eta)[2]+1
     alpha_hat<-varb_hat[startpos:dim(varb_hat)[1],startpos:dim(varb_hat)[2]]
     Tes=as.numeric(t(gamma_hat)%*%ginv(alpha_hat)%*%gamma_hat)

     rk=qr(alpha_hat)$rank
     rlt=pchisq(Tes,rk,lower.tail=FALSE)
     rlt
}
