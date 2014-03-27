innerEpi <-
function(pheno,gstnd,geno_expn,gname,gchr)
{
    gNumchr=length(gstnd)-1
    ph_idx=is.na(pheno)
    pheno=pheno[!ph_idx]
    pair_num=gNumchr*(gNumchr-1)/2
    oMat<-matrix(-1.0,pair_num,3)
    for (igIdx_i in 1:gNumchr)
    {
        print(paste(igIdx_i," of ",gNumchr," with other genes both on ",gchr[1]," chromsome!",sep=""))
        x_A=as.matrix(geno_expn[,(gstnd[igIdx_i]+1):gstnd[igIdx_i+1]])
        x_A=as.matrix(x_A[!ph_idx,])
        igIdx_j <-igIdx_i+1
        while (igIdx_j <= gNumchr)
        {
            x_B=as.matrix(geno_expn[,(gstnd[igIdx_j]+1):gstnd[igIdx_j+1]])
            x_B=as.matrix(x_B[!ph_idx,])
            rOP=(igIdx_i*(2*gNumchr-igIdx_i-1))/2+igIdx_j-gNumchr
            oMat[rOP,1]=gname[igIdx_i]
            oMat[rOP,2]=gname[igIdx_j]
            oMat[rOP,3]=as.numeric(try(fRGInteraction(pheno,x_A,x_B)))
            igIdx_j <-igIdx_j+1
        }
    }
    return(oMat)
}
