innerSnpListInteraction <-
function(pheno,snpList)
{
    rlt_out<-data.frame( )
    snp_num=dim(snpList)[2]
    snp_name<-colnames(snpList)
  
    for(i in 1:snp_num)
    {    
        j=i+1
        while(j<=snp_num)
        {            
            pval=snpPairInteraction(pheno,snpList[,i],snpList[,j])

            rOP=(i*(2*snp_num-i-1))/2+j-snp_num
    
            rlt_out[rOP,"snp1"] <-snp_name[i]
            rlt_out[rOP,"snp2"] <-snp_name[j]
            rlt_out[rOP,"pval"] <-pval

            j=j+1
        }
    }
    return(rlt_out)

}
