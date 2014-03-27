outerEpi <-
function(pheno,gstnd,gStndp,geno_expn,gname,gNamep,gchr,gChrp)
{
    
    gNumchr=length(gstnd)-1
    gene_num_postchr=length(gStndp)-1

    ph_idx=is.na(pheno)
    pheno=pheno[!ph_idx]
    
    oMat<-matrix(-1.0,gNumchr*gene_num_postchr,3)

    for (ogIdx_i  in 1:gNumchr)
    {
         print(paste(ogIdx_i," of ",gNumchr," on  ",gchr[1]," chromsome with other genes on ",gChrp[1]," chromsome(",gene_num_postchr,"genes)",sep=""))

        x_A=as.matrix(geno_expn[,(gstnd[ogIdx_i]+1):gstnd[ogIdx_i+1]])
        x_A=as.matrix(x_A[!ph_idx,])

        for(ogIdx_j in 1:gene_num_postchr)
        {                
            x_B=as.matrix(geno_expn[,(gStndp[ogIdx_j]+1):gStndp[ogIdx_j+1]])
            x_B=as.matrix(x_B[!ph_idx,])

            oPtr=(ogIdx_i-1)*gene_num_postchr+ogIdx_j
            oMat[oPtr,1]<-gname[ogIdx_i]
            oMat[oPtr,2]<-gNamep[ogIdx_j]
            oMat[oPtr,3]<-as.numeric(try(fRGInteraction(pheno,x_A,x_B)))

        }    
    } 
    
    return(oMat)

}
