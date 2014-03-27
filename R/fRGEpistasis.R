fRGEpistasis <-
function(wDir,phenoInfo,gnoFiles,mapFiles,gLst,fdr,rng)
{
    gInfo=reduceGeno(wDir,phenoInfo,gnoFiles,mapFiles,gLst,rng)
    pheno <-as.vector( phenoInfo[,2] )
    pheno <-pheno-mean(pheno)

    output<-matrix(-1,0,3)
    epi_out<-data.frame( )

    for(chr_idx in 1:gInfo$chr_num)
    {
        chrs<-unique(gInfo$gene_chr)            
        gidx=gInfo$gene_chr==chrs[chr_idx]
        gname<-gInfo$gene_name[gidx]
        gchr <-gInfo$gene_chr[gidx]
        idx_vec<-which(gInfo$gene_chr==chrs[chr_idx])
        gstnd<-c(gInfo$gene_stnd[idx_vec[1]],gInfo$gene_stnd[idx_vec+1]) 
                
        print(paste("Performing epistasis test inner chromsome",gchr[1]))
    
        oMat<-innerEpi(pheno,gstnd,gInfo$geno_expn,gname,gchr)    
        output<-rbind(output,oMat)
        rm(oMat)    
        chr_idx_post=chr_idx+1
        while (chr_idx_post <= gInfo$chr_num)
        { 
        
            gidx_post=gInfo$gene_chr==chrs[chr_idx_post]

            gNamep<-gInfo$gene_name[gidx_post]

            gChrp <-gInfo$gene_chr[gidx_post]
    
            idxVecp<-which(gInfo$gene_chr==chrs[chr_idx_post])


            gStndp<-c(gInfo$gene_stnd[idxVecp[1]],gInfo$gene_stnd[idxVecp+1])
    

            print(paste("Performing epistasis test outer ",gchr[1],":",gChrp[1],"chromsomes!"))
        
            oMat=outerEpi(pheno,gstnd,gStndp,gInfo$geno_expn,gname,gNamep,gchr,gChrp)

            output<-rbind(output,oMat)
        
            rm(oMat)
            chr_idx_post <-chr_idx_post+1
        }

    }
    if(fdr<1)
    {
        pval<-as.numeric(output[,3])
        qval<-pval*length(pval)/rank(pval)
        ## or use below
        #qval<-p.adjust(pval,method="fdr",length(pval))
        idx<-qval<=fdr
        fdrout<-output[idx,]
    }else
    {
        fdrout<-output
    }

    for(i in 1:dim(fdrout)[1])
    {
        epi_out[i,"Gene1"]=fdrout[i,1]
        epi_out[i,"Gene2"]=fdrout[i,2]
        epi_out[i,"FRG"]=fdrout[i,3]
    }

    return(epi_out)
}
