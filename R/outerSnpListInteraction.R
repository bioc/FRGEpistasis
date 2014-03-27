outerSnpListInteraction <-
function(pheno,snpList1,snpList2)
{
    rlt_out<-data.frame( )
    snp_num1=dim(snpList1)[2]
    snp_num2=dim(snpList2)[2]

    snp_name1<-colnames(snpList1)
    snp_name2<-colnames(snpList2)

    for(i in 1:snp_num1)
    {    
        for(j in 1:snp_num2)
        {            
            pval=snpPairInteraction(pheno,snpList1[,i],snpList2[,j])
            rOP=(i-1)*snp_num2+j
            rlt_out[rOP,"snp1"]<-snp_name1[i]
            rlt_out[rOP,"snp2"]<-snp_name2[j]
            rlt_out[rOP,"pval"]<-pval
        }
    }
    return(rlt_out)

}
