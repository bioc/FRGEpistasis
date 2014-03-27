reduceGeno <-
function(wDir,pheno,gnoFiles,mapFiles,gLst,rng)
{
    spl_num=dim(pheno)[1]
    gInfo=list()
    gInfo$chr_num=0
    gInfo$geno_expn=matrix(0,spl_num,1)
    gInfo$geno_expn=gInfo$geno_expn[,-1]
    gInfo$gene_name=rep(0,dim(gLst)[1])
    gInfo$gene_chr=rep(0,dim(gLst)[1])
    gInfo$gene_stnd=rep(0,dim(gLst)[1]+1)
    gCount=1
    gInfo$gNumPerChr=rep(0,dim(gnoFiles)[1])
    for(file_idx in 1:dim(gnoFiles)[1])
    {
        fdata<-read.table(paste(wDir,gnoFiles[file_idx,1],sep=""),header=TRUE)
        mapChr<-read.table(paste(wDir,mapFiles[file_idx,1],sep=""))

        gnoChr=as.matrix(fdata[,-1:-6 ] )
    
        cur_chrom=mapChr[1,1]        
        idx= gLst$Chromosome == cur_chrom
        gListChr <-gLst[idx,]

        gene.symbol<-as.vector( gListChr$Gene_Symbol)
        rep=FALSE
        chrs<-unique(gInfo$gene_chr)
        for(i in 1:length(chrs))
        {
            if(chrs[i]==cur_chrom)rep=TRUE
        }
        if(rep!=TRUE)
            gInfo$chr_num=gInfo$chr_num+1

        for (gIdx_i in 1:dim(gListChr)[1])
        {
            print(paste("Expansion gene",gIdx_i," of ",dim(gListChr)[1]," on chromsome",cur_chrom, "!",sep=""))
            x_A<-fourierExpansion(gIdx_i,gnoChr,gListChr,mapChr,rng)
            if(dim(x_A)[2]>0)
            {
                gInfo$geno_expn=cbind(gInfo$geno_expn,x_A)
                gInfo$gene_name[gCount]=gene.symbol[gIdx_i]
                gInfo$gene_chr[gCount]=cur_chrom
                gInfo$gene_stnd[gCount+1]=gInfo$gene_stnd[gCount]+dim(x_A)[2]
                gInfo$gNumPerChr[file_idx]=gInfo$gNumPerChr[file_idx]+1
                gCount=gCount+1

            }
                
        }
    }
    idx_0=gInfo$gene_name==0
    gInfo$gene_name<-gInfo$gene_name[!idx_0]
    nonempty_gene_num=length(gInfo$gene_name)
    gInfo$gene_chr <-gInfo$gene_chr[1:nonempty_gene_num]
    gInfo$gene_stnd<-gInfo$gene_stnd[1:(nonempty_gene_num+1)]
    
    return(gInfo)
}
