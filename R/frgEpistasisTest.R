frgEpistasisTest <-
function(pheno,geno_A,pos_A,geno_B,pos_B)
{    

    sample_num<-dim(geno_A)[1]
    pheno<-pheno-mean(pheno)
    
    if( !is.null( dim(geno_A) ) && dim(geno_A)[2]>3)
    {
        maf<-colMeans(geno_A,na.rm=TRUE)/2
        geno_A[,maf>0.5]=2-geno_A[,maf>0.5]
        geno_A[is.na(geno_A)]=0
        maf=colMeans(geno_A)/2
        pos_A<-pos_A[maf>0]
        geno_A<-geno_A[,maf>0]
        maf=maf[maf>0]
      }
    snp_num_A=dim(geno_A)[2]
            
    if(snp_num_A>3)
    {
        if ( is.null(pos_A) )
        {
            pos_A <-(0:( snp_num_A-1) )/(snp_num_A-1)
        }else {
            idx<-order(pos_A)
            geno_A<-geno_A[,idx]
            pos_A<-pos_A[idx]
            pos_A<-(pos_A-pos_A[1])/(pos_A[snp_num_A]-pos_A[1])
        }

        # use PCA to determine the number of basis
        eigenval<-prcomp(geno_A)$sd^2
        sum_eigen=sum(eigenval)
        tmp=0
        n_of_basis=0
        for(i in 1:length(eigenval))
        {
            tmp=eigenval[i]+tmp
            n_of_basis=i;
            if(tmp>=0.8*sum_eigen) break
        }
        n_of_basis=floor(n_of_basis/2)*2+1 
        #make n_of_basis_A the odd number
        # end of setting the basis number

        frange <-c(pos_A[1], pos_A[length(pos_A)])
        
        fbasis<-create.fourier.basis(frange,nbasis=n_of_basis)
        phi=eval.basis(pos_A,fbasis);
        x_A <-t(ginv(t(phi)%*%phi)%*%t(phi)%*%t(geno_A))
    }else x_A  <-geno_A

    
    if( !is.null( dim(geno_B) ) && dim(geno_B)[2]>3)
    {
        maf<-colMeans(geno_B,na.rm=TRUE)/2
        geno_B[,maf>0.5]=2-geno_B[,maf>0.5]
        geno_B[is.na(geno_B)]=0
        maf=colMeans(geno_B)/2
        pos_B<-pos_B[maf>0]
        geno_B<-geno_B[,maf>0]
        maf=maf[maf>0]
      }
    snp_num_B=dim(geno_B)[2]
    
    if(snp_num_B>3)
    {
        if ( is.null(pos_B) )
        {
            pos_B <-(0:( snp_num_B-1) )/(snp_num_B-1)
        }else {
            idx<-order(pos_B)
            geno_B<-geno_B[,idx]
            pos_B<-pos_B[idx]
            pos_B<-(pos_B-pos_B[1])/(pos_B[snp_num_B]-pos_B[1])
        }

        # use PCA to determine the number of basis
        eigenval<-prcomp(geno_B)$sd^2
        sum_eigen=sum(eigenval)
        tmp=0
        n_of_basis=0
        for(i in 1:length(eigenval))
        {
            tmp=eigenval[i]+tmp
            n_of_basis=i;
            if(tmp>=0.8*sum_eigen) break
        }
        n_of_basis=floor(n_of_basis/2)*2+1 
        #make n_of_basis_A the odd number
        # end of setting the basis number

        frange <-c(pos_B[1], pos_B[length(pos_B)])
        
        fbasis<-create.fourier.basis(frange,nbasis=n_of_basis)
        phi=eval.basis(pos_B,fbasis);
        x_B <-t(ginv(t(phi)%*%phi)%*%t(phi)%*%t(geno_B))
    }else x_B<-geno_B
    
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
     b_hat=iWTW%*%t(W)%*%pheno
      
     theta<-1/(length(pheno)-basis_num_A-basis_num_B-basis_num_A*basis_num_B)
     delta=theta*as.numeric(t(pheno)%*%pheno-t(pheno)%*%W%*%b_hat)
     varb_hat=delta*iWTW

     gamma_hat<-b_hat[(basis_num_A+basis_num_B+1):length(b_hat),]
     startpos<-basis_num_A+basis_num_B+1
     alpha_hat<-varb_hat[startpos:dim(varb_hat)[1],startpos:dim(varb_hat)[2]]
     Tes=as.numeric(t(gamma_hat)%*%ginv(alpha_hat)%*%gamma_hat)

     rk=qr(alpha_hat)$rank
     rlt=pchisq(Tes,rk,lower.tail=FALSE)

     rlt

}
