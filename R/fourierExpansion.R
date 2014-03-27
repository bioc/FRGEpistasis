fourierExpansion <-
function(gene_idx,geno,gene_list,snp_map,rng)
{
    
    snpName<-as.vector( snp_map[,2] )
    snp.position<-as.vector( snp_map[,4] )
    snp.chrom=snp_map[,1]
    
    gene.chrom<-as.vector( gene_list$Chromosome )
    gene.symbol<-as.vector( gene_list$Gene_Symbol)
    gene.start<-as.vector(gene_list$Start)
    gene.end<-as.vector(gene_list$End)

    sub.idx<-snp.chrom == gene.chrom[gene_idx]
    sub.idx<-sub.idx & (snp.position>(gene.start[gene_idx] -rng) )
    sub.idx<-sub.idx & (snp.position<(gene.end[gene_idx]+rng) )
    x<-as.matrix(geno[,sub.idx])
    pos<-snp.position[sub.idx]
    
    if( !is.null( dim(x) ))
    {
        maf<-colMeans(x,na.rm=TRUE)/2
        x[,maf>0.5]=2-x[,maf>0.5]
        x[is.na(x)]=0
        maf=colMeans(x)/2
        pos<-pos[maf>0]
        pos=(pos-pos[1])/(pos[length(pos)]-pos[1])
        x<-as.matrix(x[,maf>0])
        maf=maf[maf>0]
    }
        nsnps <-dim(x)[2]
    
    gil=list();
    # if the number of snp in the gene is over 3, expansion occurs
    if(nsnps>3)
    {
        if ( is.null(pos) )
        {
            pos <-(0:( nsnps-1) )/(nsnps-1)
        }else {
            idx<-order(pos)
            x<-x[,idx]
            pos<-pos[idx]
            pos<-(pos-pos[1])/(pos[nsnps]-pos[1])
        }

        # use PCA to determine the number of basis
        eigenval<-prcomp(x)$sd^2
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

        frange <-c(pos[1], pos[length(pos)])
        
        fbasis<-create.fourier.basis(frange,nbasis=n_of_basis)
        phi=eval.basis(pos,fbasis);
        gil$zeta <-t(ginv(t(phi)%*%phi)%*%t(phi)%*%t(x))
    }else gil$zeta  <-x

    return(gil$zeta)

}
