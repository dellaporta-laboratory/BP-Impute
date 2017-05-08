args <- commandArgs(trailingOnly = TRUE)

vcf_path <- args[1]
skips <- args[2]
IR64 <- args[3]
DivDonor <- args[4]
hetparstat <- as.numeric(args[5])
trimstat <- args[6]
keephet <- args[7]


print('Processing vcf file and marker input.')

#if there is an error, it might be bc the number of skips is wrong
nam=read.table(vcf_path,header = T,skip = skips,as.is=T,comment="")
colnames(nam)[1]="CHROM"

nam_name=unlist(strsplit(unlist(strsplit(vcf_path,'/'))[6],'_'))[1]

#10: because we just want the genotypes, no metadata
genos_only=apply(X = nam[,10:ncol(nam)],MARGIN = c(1,2),function(x){unlist(strsplit(x,':'))[1]})

#standardize heterozygosity, just in case
genos_only[genos_only=='1/0']='0/1'

#keep genotypes for emission probabilities
genos_original=genos_only

#convert genotypes to dosage
IR64_match_mat=apply(genos_only,2,function(x){x==genos_only[,IR64]})
het_match_mat=apply(genos_only,2,function(x){x==rep('0/1',nrow(genos_only))})
DivDonor_match_mat=apply(genos_only,2,function(x){x==genos_only[,DivDonor]})

table(genos_only)

#DivDon dosage
genos_only[IR64_match_mat]=0
genos_only[het_match_mat]=1
genos_only[DivDonor_match_mat]=2
genos_only[genos_only=='./.']=NA

dim(genos_only)
dim(nam)
table(genos_only)

genos_only=matrix(as.numeric(genos_only),nrow=nrow(genos_only))
colnames(genos_only)=colnames(nam[,10:ncol(nam)])

#if hetparstat is true, remove parentals and highly heterozygous lines
if(hetparstat>0){
  het_prop=unlist(apply(genos_only,2,function(x) sum(x==1,na.rm = T)/sum(!is.na(x))))
  hetparkeep=het_prop<hetparstat & !(colnames(genos_only) %in% c(IR64,DivDonor))
  names(hetparkeep)=colnames(genos_only)
  
  #output lines that are kept for archival purposes
  write.table(t(hetparkeep),paste(nam_name,'_',hetparstat,'_het_pars_removed.txt',sep=''),col.names = T,quote=F,sep = '\t')
  
  #remove lines from data
  #for the nam file, we also need to keep the first 9 metadata columns
  nam=nam[,c(rep(TRUE,9),hetparkeep)]
  genos_only=genos_only[,hetparkeep]
}

#Now, if we are to keep any heterozygosity at all (for the genetic map)
#Any heterozygosity will be imputed as missing
genos_opts=c(0,0.5,1)
if(keephet=='nohet'){
  nam[,10:ncol(nam)][genos_only==1]='./.'
  genos_only[genos_only==1]=NA  
  genos_opts=c(0,1)
}

#we also need a matrix of coverage
#coverage=apply(X = nam[,10:ncol(nam)],MARGIN = c(1,2),function(x){as.numeric(unlist(strsplit(unlist(strsplit(x,':'))[2],',')))})

ref_coverage=apply(X = nam[,10:ncol(nam)],MARGIN = c(1,2),function(x){unlist(strsplit(unlist(strsplit(x,':'))[2],','))[1]})
alt_coverage=apply(X = nam[,10:ncol(nam)],MARGIN = c(1,2),function(x){unlist(strsplit(unlist(strsplit(x,':'))[2],','))[2]})

#There will be coerced NAs in ref_coverage due to './.'
ref_coverage=matrix(as.numeric(ref_coverage),nrow=nrow(ref_coverage))
alt_coverage=matrix(as.numeric(alt_coverage),nrow=nrow(alt_coverage))

#set NA's in coverage to 0
ref_coverage[is.na(ref_coverage)]=0
alt_coverage[is.na(alt_coverage)]=0

#imputed marker indicator matrices
#NA's from missing markers were already set above to be have a coverage of 0
sequed_mat=ref_coverage!=99 & alt_coverage!=99 & (ref_coverage!=0 | alt_coverage!=0)
imp_mat=ref_coverage==99 & alt_coverage==99
impseq_mat=sequed_mat | imp_mat

#forward_pass_true_seqimp=genos_only
#reverse_pass_true_seqimp=genos_only


#remove parentals
#full_genos_only[,!(colnames(full_genos_only) %in% IR64)]=NULL
#full_genos_only[,!(colnames(full_genos_only) %in% DivDonor)]=NULL


#set imputed=1, missing probabilities = NA
#update of sequencing marker binomial probabilities will be calculated given the parental state of forward pass or the actual parental state of the marker
forward_probs=imp_mat
reverse_probs=imp_mat

#matrices that represent imputation if the forward pass were all correct, or if the reverse pass were all correct
forward_pass_true=genos_only
forward_pass_true[sequed_mat]=NA

reverse_pass_true=genos_only
reverse_pass_true[sequed_mat]=NA

chroms=names(table(nam$CHROM))
#chroms=c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10','Chr11','Chr12')

print('Imputing...')
for(chrom in chroms){
  temp_forward_pass_true=forward_pass_true[nam$CHROM==chrom,]
  temp_reverse_pass_true=reverse_pass_true[nam$CHROM==chrom,]
  temp_genos_only=genos_only[nam$CHROM==chrom,]
  temp_forward_probs=forward_probs[nam$CHROM==chrom,]
  temp_reverse_probs=reverse_probs[nam$CHROM==chrom,]
  temp_sequed_mat=sequed_mat[nam$CHROM==chrom,]
  temp_alt_coverage=alt_coverage[nam$CHROM==chrom,]
  temp_ref_coverage=ref_coverage[nam$CHROM==chrom,]
  temp_genos_original=genos_original[nam$CHROM==chrom,]
  temp_imp_mat=imp_mat[nam$CHROM==chrom,]

  for(j in 1:ncol(temp_genos_only)){
  #for(j in 107:107){
    print(c(chrom,j,nam_name))
    #forward pass
    i=2
    while(i <= nrow(temp_forward_pass_true)){
      #print(c(j,i))
      #identify missing interval
      if(is.na(temp_forward_pass_true[i,j])){
        #print('yee')
        tm1=i
        t=i
        while(is.na(temp_forward_pass_true[tm1,j])){if(tm1==1){break};tm1=tm1-1;}
        while(is.na(temp_forward_pass_true[t,j])){if(t==nrow(temp_forward_pass_true)){break};t=t+1;} 
        
        #what is the best way to measure this?
        #remember, there is very little heterozygosity, especially after removing heterozygous lines
        p_total=sum(temp_genos_only[tm1,]!=temp_genos_only[t,],na.rm = T)/sum(!is.na(temp_genos_only[tm1,]) & !is.na(temp_genos_only[t,]))
        #print(p_total)
        if(p_total==0){p_total=1/1000}
        
        #incorporate new probabilities and genotypes in missing interval
        x=tm1+1
        while(x<=t){

          #impute sequenced markers probabilities (not genotypes); emission probabilities conditional on previous parental states
          #what is the probability of this marker if we were going to impute it from the previous parental state
          #push the previous chain up one more marker to get a mixture of probabilities at sequenced markers
          emission=1
                    
          #in case the first marker is sequed, set the prior probability here
          if((i-1)==1 & temp_sequed_mat[i-1,j]){
            if(temp_genos_only[i-1,j]==0){
              if(temp_genos_original[i-1,IR64]=='0/0'){emission=(1-0.05)^temp_ref_coverage[i-1,j]*0.05^(temp_alt_coverage[i-1,j])}
              if(temp_genos_original[i-1,IR64]=='1/1'){emission=(1-0.05)^temp_alt_coverage[i-1,j]*0.05^(temp_ref_coverage[i-1,j])}
            }
            if(temp_genos_only[i-1,j]==2){
              if(temp_genos_original[i-1,DivDonor]=='0/0'){emission=(1-0.05)^temp_ref_coverage[i-1,j]*0.05^(temp_alt_coverage[i-1,j])}
              if(temp_genos_original[i-1,DivDonor]=='1/1'){emission=(1-0.05)^temp_alt_coverage[i-1,j]*0.05^(temp_ref_coverage[i-1,j])}
            }
            if(temp_genos_only[i-1,j]==1){
              emission=(0.5)^temp_ref_coverage[i-1,j]*0.5^(temp_alt_coverage[i-1,j])
            } 
            #use the sequenced marker as a distal anchor
            temp_forward_probs[i-1,j]=emission
            #print(c(as.numeric(temp_forward_probs[i-1-1,j]),(1-p_total),emission))
            #impute the missing and sequenced markers
            temp_forward_pass_true[i-1,j]=temp_genos_only[i-1,j]          
          }
          
          if(temp_imp_mat[x,j]){
            temp_forward_probs[x,j]=1
            #print(c(as.numeric(temp_forward_probs[x-1,j]),(1-p_total),emission))
            #impute the missing and sequenced markers
            temp_forward_pass_true[x,j]=temp_genos_only[x,j]
            x=x+1
            next
          }
          
          if(as.numeric(temp_forward_probs[x-1,j])<1e-99 & as.numeric(temp_forward_probs[x-1,j])>0){
            temp_forward_probs[x,j]=1e-100
            temp_forward_pass_true[x,j]=temp_forward_pass_true[x-1,j]
            x=x+1
            next
          } 
          
          if(!temp_sequed_mat[x,j]){
            temp_forward_probs[x,j]=as.numeric(temp_forward_probs[x-1,j])*(1-p_total)*emission
            #print(c(as.numeric(temp_forward_probs[x-1,j]),(1-p_total),emission))
            #impute the missing and sequenced markers
            temp_forward_pass_true[x,j]=temp_forward_pass_true[x-1,j]
            x=x+1
            next
          }
          if(!is.na(temp_forward_pass_true[x-1,j]) & temp_sequed_mat[x,j]){
            if(temp_forward_pass_true[x-1,j]==0){
              if(temp_genos_original[x,IR64]=='0/0'){emission=(1-0.05)^temp_ref_coverage[x,j]*0.05^(temp_alt_coverage[x,j])}
              if(temp_genos_original[x,IR64]=='1/1'){emission=(1-0.05)^temp_alt_coverage[x,j]*0.05^(temp_ref_coverage[x,j])}
            }
            if(temp_forward_pass_true[x-1,j]==2){
              if(temp_genos_original[x,DivDonor]=='0/0'){emission=(1-0.05)^temp_ref_coverage[x,j]*0.05^(temp_alt_coverage[x,j])}
              if(temp_genos_original[x,DivDonor]=='1/1'){emission=(1-0.05)^temp_alt_coverage[x,j]*0.05^(temp_ref_coverage[x,j])}
            }
            if(temp_forward_pass_true[x-1,j]==1){
              emission=(0.5)^temp_ref_coverage[x,j]*0.5^(temp_alt_coverage[x,j])
            }
            temp_forward_probs[x,j]=as.numeric(temp_forward_probs[x-1,j])*(1-p_total)*emission
            #print(c(as.numeric(temp_forward_probs[x-1,j]),(1-p_total),emission))
            #impute the missing and sequenced markers
            temp_forward_pass_true[x,j]=temp_forward_pass_true[x-1,j]
            x=x+1
            next
          }
          
          if(is.na(temp_forward_pass_true[x-1,j]) & temp_sequed_mat[x,j]){
            if(temp_genos_only[x,j]==0){
              if(temp_genos_original[x,IR64]=='0/0'){emission=(1-0.05)^temp_ref_coverage[x,j]*0.05^(temp_alt_coverage[x,j])}
              if(temp_genos_original[x,IR64]=='1/1'){emission=(1-0.05)^temp_alt_coverage[x,j]*0.05^(temp_ref_coverage[x,j])}
            }
            if(temp_genos_only[x,j]==2){
              if(temp_genos_original[x,DivDonor]=='0/0'){emission=(1-0.05)^temp_ref_coverage[x,j]*0.05^(temp_alt_coverage[x,j])}
              if(temp_genos_original[x,DivDonor]=='1/1'){emission=(1-0.05)^temp_alt_coverage[x,j]*0.05^(temp_ref_coverage[x,j])}
            }
            if(temp_genos_only[x,j]==1){
              emission=(0.5)^temp_ref_coverage[x,j]*0.5^(temp_alt_coverage[x,j])
            } 
            #use the sequenced marker as a distal anchor
            temp_forward_probs[x,j]=emission
            #print(c(as.numeric(temp_forward_probs[x-1,j]),(1-p_total),emission))
            #impute the missing and sequenced markers
            temp_forward_pass_true[x,j]=temp_genos_only[x,j] 
            x=x+1
          }          
        }
        i=t
      }
      i=i+1  
    }
    #reverse pass
   i=(nrow(temp_reverse_pass_true)-1)
   while(i >= 1){
      #print(c(j,i))
      if(is.na(temp_reverse_pass_true[i,j])){
        tm1=i
        t=i
        while(is.na(temp_reverse_pass_true[tm1,j])){if(tm1==1){break};tm1=tm1-1;}
        while(is.na(temp_reverse_pass_true[t,j])){if(t==nrow(temp_reverse_pass_true)){break};t=t+1;} 
        
        #what is the best way to measure this?
        p_total=sum(temp_genos_only[tm1,]!=temp_genos_only[t,],na.rm = T)/sum(!is.na(temp_genos_only[tm1,]) & !is.na(temp_genos_only[t,]))
        #print(p_total)
        if(p_total==0){p_total=1/1000}
        
        x=t-1
        while(x>=tm1){ 

          #impute sequenced markers probabilities (not genotypes); emission probabilities conditional on previous parental states
          #what is the probability of this marker if we were going to impute it from the previous parental state
          #push the previous chain up one more marker to get a mixture of probabilities at sequenced markers
          emission=1
          
          #in case the first marker is sequed, set the prior probability here
          if((i+1)==nrow(temp_reverse_pass_true) & temp_sequed_mat[i+1,j]){
            if(temp_genos_only[i+1,j]==0){
              if(temp_genos_original[i+1,IR64]=='0/0'){emission=(1-0.05)^temp_ref_coverage[i+1,j]*0.05^(temp_alt_coverage[i+1,j])}
              if(temp_genos_original[i+1,IR64]=='1/1'){emission=(1-0.05)^temp_alt_coverage[i+1,j]*0.05^(temp_ref_coverage[i+1,j])}
            }
            if(temp_genos_only[i+1,j]==2){
              if(temp_genos_original[i+1,DivDonor]=='0/0'){emission=(1-0.05)^temp_ref_coverage[i+1,j]*0.05^(temp_alt_coverage[i+1,j])}
              if(temp_genos_original[i+1,DivDonor]=='1/1'){emission=(1-0.05)^temp_alt_coverage[i+1,j]*0.05^(temp_ref_coverage[i+1,j])}
            }
            if(temp_genos_only[i+1,j]==1){
              emission=(0.5)^temp_ref_coverage[i+1,j]*0.5^(temp_alt_coverage[i+1,j])
            } 
            #use the sequenced marker as a distal anchor
            temp_reverse_probs[i+1,j]=emission
            #print(c(as.numeric(temp_reverse_probs[i+1,j]),(1-p_total),emission))
            #impute the missing and sequenced markers
            temp_reverse_pass_true[i+1,j]=temp_genos_only[i+1,j]         
          }
          
          if(temp_imp_mat[x,j]){
            temp_reverse_probs[x,j]=1
            #print(c(as.numeric(temp_forward_probs[x-1,j]),(1-p_total),emission))
            #impute the missing and sequenced markers
            temp_reverse_pass_true[x,j]=temp_genos_only[x,j]
            x=x-1
            next
          }
          
          if(as.numeric(temp_reverse_probs[x+1,j])<1e-99 & as.numeric(temp_reverse_probs[x+1,j])>0){
            temp_reverse_probs[x,j]=1e-100
            temp_reverse_pass_true[x,j]=temp_reverse_pass_true[x+1,j]
            x=x-1
            next
          }  
          
          if(!temp_sequed_mat[x,j]){
            temp_reverse_probs[x,j]=as.numeric(temp_reverse_probs[x+1,j])*(1-p_total)*emission
            #print(c(as.numeric(temp_reverse_probs[x-1,j]),(1-p_total),emission))
            #impute the missing and sequenced markers
            temp_reverse_pass_true[x,j]=temp_reverse_pass_true[x+1,j]
            x=x-1
            next
          }
          if(!is.na(temp_reverse_pass_true[x+1,j]) & temp_sequed_mat[x,j]){
            if(temp_reverse_pass_true[x+1,j]==0){
              if(temp_genos_original[x,IR64]=='0/0'){emission=(1-0.05)^temp_ref_coverage[x,j]*0.05^(temp_alt_coverage[x,j])}
              if(temp_genos_original[x,IR64]=='1/1'){emission=(1-0.05)^temp_alt_coverage[x,j]*0.05^(temp_ref_coverage[x,j])}
            }
            if(temp_reverse_pass_true[x+1,j]==2){
              if(temp_genos_original[x,DivDonor]=='0/0'){emission=(1-0.05)^temp_ref_coverage[x,j]*0.05^(temp_alt_coverage[x,j])}
              if(temp_genos_original[x,DivDonor]=='1/1'){emission=(1-0.05)^temp_alt_coverage[x,j]*0.05^(temp_ref_coverage[x,j])}
            }
            if(temp_reverse_pass_true[x+1,j]==1){
              emission=(0.5)^temp_ref_coverage[x,j]*0.5^(temp_alt_coverage[x,j])
            }
            temp_reverse_probs[x,j]=as.numeric(temp_reverse_probs[x+1,j])*(1-p_total)*emission
            #print(c(as.numeric(temp_reverse_probs[x-1,j]),(1-p_total),emission))
            #impute the missing and sequenced markers
            temp_reverse_pass_true[x,j]=temp_reverse_pass_true[x+1,j]
            x=x-1
            next
          }
          
          if(is.na(temp_reverse_pass_true[x+1,j]) & temp_sequed_mat[x,j]){
            if(temp_genos_only[x,j]==0){
              if(temp_genos_original[x,IR64]=='0/0'){emission=(1-0.05)^temp_ref_coverage[x,j]*0.05^(temp_alt_coverage[x,j])}
              if(temp_genos_original[x,IR64]=='1/1'){emission=(1-0.05)^temp_alt_coverage[x,j]*0.05^(temp_ref_coverage[x,j])}
            }
            if(temp_genos_only[x,j]==2){
              if(temp_genos_original[x,DivDonor]=='0/0'){emission=(1-0.05)^temp_ref_coverage[x,j]*0.05^(temp_alt_coverage[x,j])}
              if(temp_genos_original[x,DivDonor]=='1/1'){emission=(1-0.05)^temp_alt_coverage[x,j]*0.05^(temp_ref_coverage[x,j])}
            }
            if(temp_genos_only[x,j]==1){
              emission=(0.5)^temp_ref_coverage[x,j]*0.5^(temp_alt_coverage[x,j])
            } 
            #use the sequenced marker as a distal anchor
            temp_reverse_probs[x,j]=emission
            #print(c(as.numeric(temp_reverse_probs[x-1,j]),(1-p_total),emission))
            #impute the missing and sequenced markers
            temp_reverse_pass_true[x,j]=temp_genos_only[x,j] 
            x=x-1
          }          
        }
        i=tm1
      }
      i=i-1  
    }
  }
  forward_pass_true[nam$CHROM==chrom,]=temp_forward_pass_true
  reverse_pass_true[nam$CHROM==chrom,]=temp_reverse_pass_true
  genos_only[nam$CHROM==chrom,]=temp_genos_only
  forward_probs[nam$CHROM==chrom,]=temp_forward_probs
  reverse_probs[nam$CHROM==chrom,]=temp_reverse_probs
}

print('Writing output...')

#probabilities and genotypes of forward and reverse passes, but normalize them to sum to 1 first
#replace any NA's with probability of 0 (such as missing, distal chromosome calls)
forward_pass_true[is.na(forward_pass_true)]=0
forward_probs[is.na(forward_probs)]=0
reverse_pass_true[is.na(reverse_pass_true)]=0
reverse_probs[is.na(reverse_probs)]=0

#now sum the probabilities of 0,1,2 genotypes
zeros_matrix=matrix(rep(0,length=length(genos_only)),nrow=nrow(genos_only),ncol=ncol(genos_only))
ones_matrix=matrix(rep(0,length=length(genos_only)),nrow=nrow(genos_only),ncol=ncol(genos_only))
twos_matrix=matrix(rep(0,length=length(genos_only)),nrow=nrow(genos_only),ncol=ncol(genos_only))


zeros_matrix[forward_pass_true==0]=zeros_matrix[forward_pass_true==0]+forward_probs[forward_pass_true==0]
zeros_matrix[reverse_pass_true==0]=zeros_matrix[reverse_pass_true==0]+reverse_probs[reverse_pass_true==0]

ones_matrix[forward_pass_true==1]=ones_matrix[forward_pass_true==1]+forward_probs[forward_pass_true==1]
ones_matrix[reverse_pass_true==1]=ones_matrix[reverse_pass_true==1]+reverse_probs[reverse_pass_true==1]

twos_matrix[forward_pass_true==2]=twos_matrix[forward_pass_true==2]+forward_probs[forward_pass_true==2]
twos_matrix[reverse_pass_true==2]=twos_matrix[reverse_pass_true==2]+reverse_probs[reverse_pass_true==2]


sum_probs=(zeros_matrix+ones_matrix+twos_matrix)

#outputs
#0 to 1 diversity donor dosage
#0 to 2 diversity donor dosage
#binary matrix of imputed genotypes

# #first, cut off distal ends based on most proximal LB-Impute marker of all samples
# #forward
rownames(impseq_mat)=1:nrow(impseq_mat)
impseq_mat_rownames_keep=rownames(impseq_mat)

if(trimstat=='trim'){
  for(chrom in chroms){
    imp_forwardtest=rep(0,ncol(impseq_mat[nam$CHROM==chrom,]))
    temp_impseq_mat=impseq_mat[nam$CHROM==chrom,]
    for(i in rownames(temp_impseq_mat)){
      imp_forwardtest=imp_forwardtest+temp_impseq_mat[i,]
      if((0 %in% imp_forwardtest)==F){break}
      #if((0 %in% temp_impseq_mat[i,])==F){break}
      impseq_mat_rownames_keep=impseq_mat_rownames_keep[impseq_mat_rownames_keep!=i]
    }
  }
  
  #reverse
  for(chrom in chroms){
    imp_reversetest=rep(0,ncol(impseq_mat[nam$CHROM==chrom,]))
    temp_impseq_mat=impseq_mat[nam$CHROM==chrom,]
    for(i in rev(rownames(temp_impseq_mat))){
      imp_reversetest=imp_reversetest+temp_impseq_mat[i,]
      if((0 %in% imp_reversetest)==F){break}
      #if((0 %in% temp_impseq_mat[i,])==F){break}
      impseq_mat_rownames_keep=impseq_mat_rownames_keep[impseq_mat_rownames_keep!=i]
    }
  }  
}
#normalize the probabilities to sum to 1 
zeros_matrix=zeros_matrix/sum_probs
ones_matrix=ones_matrix/sum_probs
twos_matrix=twos_matrix/sum_probs

#calculate the weighted genotypes
weighted_genos=zeros_matrix*0+ones_matrix*1+twos_matrix*2

#divide by 2 so that all genos are between 0 and 1 (they could be understood as probabilities as well)
weighted_genos=weighted_genos/2

#cut off the distal ends of the chromosomes, as calculated above
#for both genotype output and imputed binary matrix
if(trimstat=='trim'){
  #switch back to numeric so the row indexing is universal
  impseq_mat_rownames_keep=as.numeric(impseq_mat_rownames_keep)
  weighted_genos_out=weighted_genos[impseq_mat_rownames_keep,]
  colnames(weighted_genos_out)=colnames(nam)[10:ncol(nam)]
  imp_mat_out=imp_mat[impseq_mat_rownames_keep,]
  colnames(imp_mat_out)=colnames(nam)[10:ncol(nam)]
  nam=nam[impseq_mat_rownames_keep,]
  genos_only=genos_only[impseq_mat_rownames_keep,]
  forward_pass_true=forward_pass_true[impseq_mat_rownames_keep,]
  reverse_pass_true=reverse_pass_true[impseq_mat_rownames_keep,]
  forward_probs=forward_probs[impseq_mat_rownames_keep,]
  reverse_probs=reverse_probs[impseq_mat_rownames_keep,]
  
  #assign genotypes and imputed status to the distal markers
  for(chrom in chroms){
    temp_weighted_genos_out=weighted_genos_out[nam$CHROM==chrom,]
    temp_imp_mat_out=imp_mat_out[nam$CHROM==chrom,]
    temp_forward_pass_true=forward_pass_true[nam$CHROM==chrom,]
    temp_reverse_pass_true=reverse_pass_true[nam$CHROM==chrom,]
    temp_genos_only=genos_only[nam$CHROM==chrom,]
    dim(temp_weighted_genos_out)
    print(temp_weighted_genos_out[1,])
    print(temp_weighted_genos_out[nrow(temp_weighted_genos_out),])
    
    #set distal ends as anchors for discrete genotype assignment
    #replace the nan's (1 missing marker at distal end)
    temp_weighted_genos_out[1,][is.nan(temp_weighted_genos_out[1,])]=
      temp_weighted_genos_out[2,][is.nan(temp_weighted_genos_out[1,])]
    temp_imp_mat_out[1,]=rep(1,ncol(temp_imp_mat_out))
    
    last=nrow(temp_weighted_genos_out)
    temp_weighted_genos_out[last,][is.nan(temp_weighted_genos_out[last,])]=
      temp_weighted_genos_out[last-1,][is.nan(temp_weighted_genos_out[last,])]
    temp_imp_mat_out[last,]=rep(1,ncol(temp_imp_mat_out))
    
    #if there is a sequenced marker at the distal end of the chrom,
    #evaluate whether the forward or reverse state should be chosen
    #row 1
    for(j in 1:ncol(temp_weighted_genos_out)){
      print(c(chrom,j))
      if(temp_weighted_genos_out[1,j] %in% genos_opts){next}
      state1=temp_forward_pass_true[1,j]/2
      state2=temp_reverse_pass_true[1,j]/2
      temp_imp_mat_out[1,j]=1
            
      if((temp_weighted_genos_out[1,j]-state2)^2<(temp_weighted_genos_out[1,j]-state1)^2){temp_weighted_genos_out[1,j]=state2}
      if((temp_weighted_genos_out[1,j]-state2)^2>(temp_weighted_genos_out[1,j]-state1)^2){temp_weighted_genos_out[1,j]=state1}
      if((temp_weighted_genos_out[1,j]-state2)^2==(temp_weighted_genos_out[1,j]-state1)^2){temp_weighted_genos_out[1,j]=sample(c(state1,state2),1)} 
    }
    
    #last row
    for(j in 1:ncol(temp_weighted_genos_out)){
      last=nrow(temp_weighted_genos_out)
      if(temp_weighted_genos_out[last,j] %in% genos_opts){next}
      state1=temp_reverse_pass_true[last,j]/2
      state2=temp_forward_pass_true[last,j]/2
      temp_imp_mat_out[last,j]=1
      
      if((temp_weighted_genos_out[last,j]-state2)^2<(temp_weighted_genos_out[last,j]-state1)^2){temp_weighted_genos_out[last,j]=state2}
      if((temp_weighted_genos_out[last,j]-state2)^2>(temp_weighted_genos_out[last,j]-state1)^2){temp_weighted_genos_out[last,j]=state1}
      if((temp_weighted_genos_out[last,j]-state2)^2==(temp_weighted_genos_out[last,j]-state1)^2){temp_weighted_genos_out[last,j]=sample(c(state1,state2),1)} 
    }
    
    weighted_genos_out[nam$CHROM==chrom,]=temp_weighted_genos_out
    imp_mat_out[nam$CHROM==chrom,]=temp_imp_mat_out    
    
    #If the distal marker is a probability, apply
    #least squares and assign a genotype. Make this an anchor marker  
  }
    
  write.table(cbind(nam[,1:9],weighted_genos_out),paste(nam_name,'_weighted_genos_',trimstat,'_',keephet,'.txt',sep=''),row.names = F,col.names = T,quote = F,sep = '\t')
  write.table(cbind(nam[,1:9],imp_mat_out),paste(nam_name,'_imputed_binary_',trimstat,'_',keephet,'.txt',sep=''),row.names = F,col.names = T,quote = F,sep = '\t')
}

#don't trim if there's excessive missing data
#otherwise some of the smaller chromosomes will be entirely removed
if(trimstat!='trim'){
  weighted_genos_out=weighted_genos
  colnames(weighted_genos_out)=colnames(nam)[10:ncol(nam)]
  imp_mat_out=imp_mat
  colnames(imp_mat_out)=colnames(nam)[10:ncol(nam)]

  #assign genotypes and imputed status to the distal markers
  for(chrom in chroms){
    temp_weighted_genos_out=weighted_genos_out[nam$CHROM==chrom,]
    temp_imp_mat_out=imp_mat_out[nam$CHROM==chrom,]
    temp_forward_pass_true=forward_pass_true[nam$CHROM==chrom,]
    temp_reverse_pass_true=reverse_pass_true[nam$CHROM==chrom,]
    temp_genos_only=genos_only[nam$CHROM==chrom,]
    dim(temp_weighted_genos_out)
    print(temp_weighted_genos_out[1,])
    print(temp_weighted_genos_out[nrow(temp_weighted_genos_out),])
    
    #set distal ends as anchors for discrete genotype assignment
    #replace the nan's (1 missing marker at distal end)
    temp_weighted_genos_out[1,][is.nan(temp_weighted_genos_out[1,])]=
      temp_weighted_genos_out[2,][is.nan(temp_weighted_genos_out[1,])]
    temp_imp_mat_out[1,]=rep(1,ncol(temp_imp_mat_out))
    
    last=nrow(temp_weighted_genos_out)
    temp_weighted_genos_out[last,][is.nan(temp_weighted_genos_out[last,])]=
      temp_weighted_genos_out[last-1,][is.nan(temp_weighted_genos_out[last,])]
    temp_imp_mat_out[last,]=rep(1,ncol(temp_imp_mat_out))
    
    #if there is a sequenced marker at the distal end of the chrom,
    #evaluate whether the sequed state or imputed state should be chosen
    #row 1
    for(j in 1:ncol(temp_weighted_genos_out)){
      print(c(chrom,j))
      if(temp_weighted_genos_out[1,j] %in% genos_opts){next}
      state1=temp_forward_pass_true[1,j]/2
      state2=temp_reverse_pass_true[1,j]/2
      temp_imp_mat_out[1,j]=1
      
      if((temp_weighted_genos_out[1,j]-state2)^2<(temp_weighted_genos_out[1,j]-state1)^2){temp_weighted_genos_out[1,j]=state2}
      if((temp_weighted_genos_out[1,j]-state2)^2>(temp_weighted_genos_out[1,j]-state1)^2){temp_weighted_genos_out[1,j]=state1}
      if((temp_weighted_genos_out[1,j]-state2)^2==(temp_weighted_genos_out[1,j]-state1)^2){temp_weighted_genos_out[1,j]=sample(c(state1,state2),1)} 
    }
    
    #last row
    for(j in 1:ncol(temp_weighted_genos_out)){
      last=nrow(temp_weighted_genos_out)
      if(temp_weighted_genos_out[last,j] %in% genos_opts){next}
      state1=temp_reverse_pass_true[last,j]/2
      state2=temp_forward_pass_true[last,j]/2
      temp_imp_mat_out[last,j]=1
      
      if((temp_weighted_genos_out[last,j]-state2)^2<(temp_weighted_genos_out[last,j]-state1)^2){temp_weighted_genos_out[last,j]=state2}
      if((temp_weighted_genos_out[last,j]-state2)^2>(temp_weighted_genos_out[last,j]-state1)^2){temp_weighted_genos_out[last,j]=state1}
      if((temp_weighted_genos_out[last,j]-state2)^2==(temp_weighted_genos_out[last,j]-state1)^2){temp_weighted_genos_out[last,j]=sample(c(state1,state2),1)} 
    }
    
    weighted_genos_out[nam$CHROM==chrom,]=temp_weighted_genos_out
    imp_mat_out[nam$CHROM==chrom,]=temp_imp_mat_out    
    
    #If the distal marker is a probability, apply
    #least squares and assign a genotype. Make this an anchor marker  
  }
  
  write.table(cbind(nam[,1:9],weighted_genos_out),paste(nam_name,'_weighted_genos_',trimstat,'_',keephet,'.txt',sep=''),row.names = F,col.names = T,quote = F,sep = '\t')
  write.table(cbind(nam[,1:9],imp_mat_out),paste(nam_name,'_imputed_binary_',trimstat,'_',keephet,'.txt',sep=''),row.names = F,col.names = T,quote = F,sep = '\t')  
}

#output the number of nonmissing samples (for recombination rate measurements)
nonmissing_forward=rep(NA,nrow(genos_only))
for(i in 2:nrow(genos_only)){
  nonmissing_forward[i]=length(genos_only[i,][!is.na(genos_only[i,]) & !is.na(genos_only[i-1,])])
}

nonmissing_reverse=rep(NA,nrow(genos_only))
for(i in (nrow(genos_only)-1):1){
  nonmissing_reverse[i]=length(genos_only[i,][!is.na(genos_only[i,]) & !is.na(genos_only[i+1,])])
}

write.table(cbind(nonmissing_forward,nonmissing_reverse),paste(nam_name,'_nonmissing_counts.txt',sep = ''),quote = F,row.names = T,col.names = F)





