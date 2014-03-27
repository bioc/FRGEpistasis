pCAPiontwiseEpistasis <-
function(wDir,oEpi,phenoInfo,gnoFiles,mapFiles,gLst,rng)
{    
    oPtr=1
    pheno <-as.vector( phenoInfo[,2] )
    pheno <-pheno-mean(pheno)
    oFlag=FALSE
    if(length(oEpi) == 0) oFlag=TRUE
    for(file_idx in 1:dim(gnoFiles)[1])
    {
        fdata<-read.table(paste(wDir,gnoFiles[file_idx,1],sep=""),header=TRUE)
        mapChr<-read.table(paste(wDir,mapFiles[file_idx,1],sep=""))    
    
        gnoChr=as.matrix(fdata[,-1:-6 ] )

        cur_chrom=mapChr[1,1]        
        idx= gLst$Chromosome == cur_chrom
        gListChr <-gLst[idx,]
        gNumchr=dim(gListChr)[1]

        print(paste("Performing epistasis test inner chromsome",cur_chrom ))

        # inner analysis
        if(1<gNumchr )
             {
            for (igIdx_i in 1:gNumchr)
            {
         
                snpName<-as.vector( mapChr[,2] )
                snp.position<-as.vector( mapChr[,4] )
                snp.chrom=mapChr[,1]
    
                gene.chrom<-as.vector( gListChr$Chromosome )
                gene.symbol<-as.vector( gListChr$Gene_Symbol)
                gene.start<-as.vector(gListChr$Start)
                gene.end<-as.vector(gListChr$End)
            
                sub.idx<-snp.chrom == gene.chrom[igIdx_i]   
                sub.idx<-sub.idx & (snp.position>(gene.start[igIdx_i] -rng))
                sub.idx<-sub.idx & (snp.position<(gene.end[igIdx_i]+rng))
                x_A<-as.matrix(gnoChr[,sub.idx])
            
                if( !is.null( dim(x_A) ))
                {
                    maf<-colMeans(x_A,na.rm=TRUE)/2
                    x_A[,maf>0.5]=2-x_A[,maf>0.5]
                    x_A[is.na(x_A)]=0
                    maf=colMeans(x_A)/2        
                     x_A<-as.matrix(x_A[,maf>0])            
                    maf=maf[maf>0]

                }
                x_A<-2-x_A


                igIdx_j <-igIdx_i+1 
                    while (igIdx_j <= gNumchr)
                { 
                
                    sub.idx<-snp.chrom == gene.chrom[igIdx_j]   
                    sub.idx<-sub.idx & (snp.position>(gene.start[igIdx_j]-rng))
                    sub.idx<-sub.idx & (snp.position<(gene.end[igIdx_j]+rng))

                    x_B <-as.matrix(gnoChr[,sub.idx])
            
                    if( !is.null( dim(x_B) ))
                    {
                        maf<-colMeans(x_B,na.rm=TRUE)/2
                        x_B[,maf>0.5]=2-x_B[,maf>0.5]
                        x_B[is.na(x_B)]=0
                        maf=colMeans(x_B)/2        
                         x_B<-as.matrix(x_B[,maf>0])        
                        maf=maf[maf>0]

                    }
                    x_B<-2-x_B

                    minPval<-pointwiseInteraction(pheno,x_A,x_B)
                    pval<-pCAInteraction(pheno,x_A,x_B)
                    if(oFlag)
                    {
                        oEpi[oPtr,"Gene1"]=gene.symbol[igIdx_i]
                        oEpi[oPtr,"Gene2"]=gene.symbol[igIdx_j]
                        oEpi[oPtr,"PCA"]=pval
                        oEpi[oPtr,"Pointwise_test"]=minPval
                    }else{ 
                        oEpi[oPtr,"PCA"]=pval
                        oEpi[oPtr,"Pointwise_test"]=minPval
                    }
                    oPtr=oPtr+1

                        igIdx_j <-igIdx_j+1 
                }    
                } 

        }
    
        # outer analysis
        fIdxp=file_idx +1
        while(fIdxp<=dim(gnoFiles)[1])
        {

            fdata<-read.table(paste(wDir,gnoFiles[fIdxp,1],sep=""),header=TRUE)
            map_post<-read.table(paste(wDir,mapFiles[fIdxp,1],sep=""))
    
            geno_post= as.matrix(fdata[,-1:-6 ] )

            post_chrom=map_post[1,1]        
            idx= gLst$Chromosome == post_chrom 
            gene_list_post <-gLst[idx,]
            gene_num_post=dim(gene_list_post)[1]

            print(paste("Performing epistasis test outer ",cur_chrom ,":",post_chrom ,"chromsomes!"))


            for (ogIdx_i  in 1:gNumchr)
            {        
                snpName<-as.vector( mapChr[,2] )
                snp.position<-as.vector( mapChr[,4] )
                snp.chrom=mapChr[,1]
    
                gene.chrom<-as.vector( gListChr $Chromosome )
                gene.symbol<-as.vector( gListChr $Gene_Symbol)
                gene.start<-as.vector(gListChr $Start)
                gene.end<-as.vector(gListChr $End)
            
                sub.idx<-snp.chrom == gene.chrom[ogIdx_i]
                sub.idx<-sub.idx & (snp.position>(gene.start[ogIdx_i]-rng))
                sub.idx<-sub.idx & (snp.position<(gene.end[ogIdx_i]+rng))
                x_A<-as.matrix(gnoChr[,sub.idx])
            
                if( !is.null( dim(x_A) ))
                {
                    maf<-colMeans(x_A,na.rm=TRUE)/2
                    x_A[,maf>0.5]=2-x_A[,maf>0.5]
                    x_A[is.na(x_A)]=0
                    maf=colMeans(x_A)/2        
                     x_A<-as.matrix(x_A[,maf>0])            
                    maf=maf[maf>0]

                }
                x_A<-2-x_A

                for(ogIdx_j in 1:gene_num_post)
                {                
                    snpName<-as.vector( map_post[,2] )
                    snp.position<-as.vector( map_post[,4] )
                    snp.chrom=map_post[,1]
    
                    gene.chrom<-as.vector( gene_list_post$Chromosome )
                    gene.symbol.post<-as.vector( gene_list_post$Gene_Symbol)
                    gene.start<-as.vector(gene_list_post$Start)
                    gene.end<-as.vector(gene_list_post$End)
            
                    sub.idx<-snp.chrom == gene.chrom[ogIdx_j]   
                    sub.idx<-sub.idx & (snp.position>(gene.start[ogIdx_j]-rng))
                    sub.idx<-sub.idx & (snp.position<(gene.end[ogIdx_j]+rng)) 

                    x_B <-as.matrix(geno_post[,sub.idx])
            
                    if( !is.null(dim(x_B)))
                    {
                        maf<-colMeans(x_B,na.rm=TRUE)/2
                        x_B[,maf>0.5]=2-x_B[,maf>0.5]
                        x_B[is.na(x_B)]=0
                        maf=colMeans(x_B)/2        
                         x_B<-as.matrix(x_B[,maf>0])        
                        maf=maf[maf>0]

                    }
                    x_B<-2-x_B

                    minPval<-pointwiseInteraction(pheno,x_A,x_B)
                    pval<-pCAInteraction(pheno,x_A,x_B)
                    
                    if(oFlag)
                    {
                        oEpi[oPtr,"Gene1"]=gene.symbol[ogIdx_i]
                        oEpi[oPtr,"Gene2"]=gene.symbol.post[ogIdx_j]
                        oEpi[oPtr,"PCA"]=pval
                        oEpi[oPtr,"Pointwise_test"]=minPval
                    }else{
                        oEpi[oPtr,"PCA"]=pval
                        oEpi[oPtr,"Pointwise_test"]=minPval
                    }

                    oPtr=oPtr+1

                }    
                } 

            fIdxp=fIdxp+1
        }
    }

    return(oEpi)
}
