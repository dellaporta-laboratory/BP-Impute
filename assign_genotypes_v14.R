args <- commandArgs(trailingOnly = TRUE)

genos_path <- args[1]
imp_path <- args[2]

nam_name=unlist(strsplit(unlist(strsplit(genos_path,'/'))[6],'_'))[1]

genos=read.table(genos_path,skip=0,header=T,as.is=T)
imputed=read.table(imp_path,skip=0,header=T,as.is=T)

#not needed in new vcf files
#set parents as automatically imputed
#IR64='ID153_H10' 
#DivDon='ID153_H11'

#imputed[,IR64]=rep(TRUE,nrow(imputed))
#imputed[,DivDon]=rep(TRUE,nrow(imputed))

#select a certain line and chromosome for testing
#imputed=imputed[genos$CHROM=='Chr6',123:125]
#genos=genos[genos$CHROM=='Chr6',123:125]
#genos_only=genos
#imputed_only=as.matrix(imputed)
#chroms='Chr6'
#assigned=genos_only
#assigned[!imputed_only]=NA
#temp_assigned=assigned
#temp_imputed=imputed_only
#temp_genos_only=genos_only
#temp_genos=genos


chroms=names(table(genos$CHROM))

#remove metadata
#imputed needs to be a matrix, not a data frame
genos_only=genos[,10:ncol(genos)]
imputed_only=as.matrix(imputed[,10:ncol(imputed)])


assigned=genos_only
#assigned[imputed_only]=genos_only[imputed_only]
assigned[!imputed_only]=NA

for(chrom in chroms){
  temp_assigned=assigned[genos$CHROM==chrom,]
  temp_imputed=imputed_only[genos$CHROM==chrom,]
  temp_genos_only=genos_only[genos$CHROM==chrom,]
  temp_genos=genos[genos$CHROM==chrom,]
  
#   states=c(0,0.5,1)
#   
#   for(j in 1:ncol(temp_genos_only)){
#     names(states)=(temp_genos_only[1,j]-states)^2
#     print(c(j,names(states),temp_genos_only[1,j],states[min(names(states))]))
#     temp_genos_only[1,j]=states[min(names(states))]
#   }
  
  
  #If the distal marker is a probability, apply
  #least squares and assign a genotype. Make this an anchor marker
  
  
  #assign missing genotypes
  for(j in 1:ncol(temp_assigned)){
  #for(j in 179:179){
    evaluate=temp_assigned[,j]
    i=1
    while(i <= nrow(temp_assigned)){
    #while(i <= 6035){ 
      if(is.na(evaluate[i])){
        print(c(chrom,j,i))
        #find flanking states
        #then recursively assign genotypes to all missing in that block
        #store probabilities, in case the
        #probabilities are all halfway between the two genotypes
        
        left=NA
        right=NA
        #now search
        for(il in i:1){left=temp_assigned[il,j];if(!is.na(left)){break}}
        for(ir in i:length(evaluate)){right=temp_assigned[ir,j];if(!is.na(right)){break}}
        
        #if the search reaches the end of the chromosome, make the left and right genotypes the same
        if(is.na(left)){left=right}
        if(is.na(right)){right=left}
              
        #clean up the entire region
        #i2=i+1
        i2=i
        while(is.na(temp_assigned[i2,j]) & i2 <= nrow(temp_assigned)){
          #end=i2
          if((temp_genos_only[i2,j]-left)^2<(temp_genos_only[i2,j]-right)^2){temp_assigned[i2,j]=left;i2=i2+1;next}
          if((temp_genos_only[i2,j]-left)^2>(temp_genos_only[i2,j]-right)^2){temp_assigned[i2,j]=right;i2=i2+1;next}
          if((temp_genos_only[i2,j]-left)^2==(temp_genos_only[i2,j]-right)^2){temp_assigned[i2,j]=sample(c(left,right),1);i2=i2+1;next}          
        }
        #print(temp_assigned[start:i2,j])
        
        #identify subinterval where there may be a region that's impossible to impute
        rows=1:nrow(temp_assigned)
        rows_interval=rows>= i & rows<=i2
        ambig=rows_interval
        
        if(left==0 & right==1 | left==1 & right==0){
          ambig=ambig & temp_genos_only[,j]==0.5         
        }
        
        if(left==0 & right==0.5 | left==0.5 & right==0){
          ambig=ambig & temp_genos_only[,j]==0.25          
        }
        
        if(left==1 & right==0.5 | left==0.5 & right==1){
          ambig=ambig & temp_genos_only[,j]==0.75          
        }
        
        if(left==right){ambig=rep(F,nrow(temp_genos_only))}
      
        #print(c(sum(ambig),i))        
        #if the region has the same genotype, choose a random breakpoint
        #print(temp_assigned[start:i2,j])
        if(sum(ambig) > 1){             
          newstart=min(rows[ambig])
          newend=max(rows[ambig])
          print(c('equal_prob_region',i,j,newstart,newend))
          #the breakpoint is assumed to be to the right of the chosen marker
          #so we don't run off the chromosome
          
          #add two sites on either end
          breakpt=sample(c('ll',newstart:newend,'rr'),1)
          #print(breakpt)

          if(breakpt=='ll'){
            temp_assigned[rownames(temp_assigned)[1:nrow(temp_assigned)>= newstart & 1:nrow(temp_assigned)<=newend],j]=right            
          }
          if(breakpt=='rr'){
            temp_assigned[rownames(temp_assigned)[1:nrow(temp_assigned)>= newstart & 1:nrow(temp_assigned)<=newend],j]=left            
          }
          if(breakpt!='ll' & breakpt!='rr'){
            breakpt=as.numeric(breakpt)
            temp_assigned[rownames(temp_assigned)[1:nrow(temp_assigned)>= newstart & 1:nrow(temp_assigned)<=breakpt],j]=left
            #we overwrite the breakpoint marker above
            #we do this so that the indexing stays in the interval
            temp_assigned[rownames(temp_assigned)[1:nrow(temp_assigned)> breakpt & 1:nrow(temp_assigned)<=newend],j]=right 
          }
        }
        #collapse NAs, treat as the same 'i'
        i=i2
        #next
      }
      i=i+1
    }
  }
  assigned[genos$CHROM==chrom,]=temp_assigned 
}  

#output
out=assigned
out[out==0]='A'
out[out==0.5]='H'
out[out==1]='B'
loci=paste(rep('*',nrow(genos)),genos[,1],rep('_',nrow(genos)),genos[,2],rep(' ',nrow(genos)),sep='')
#colnames(out)=c('marker',colnames(genos)[10:ncol(genos)])


#mapmaker format
write.table('data type ri self',paste(nam_name,'_reseq_weighted_genos_trim_yeshet_assigned_genos.txt',sep=''),
            quote=F,row.names=F,col.names=F)

write.table(t(c(ncol(genos_only),nrow(genos_only),0)),paste(nam_name,'_reseq_weighted_genos_trim_yeshet_assigned_genos.txt',sep=''),
            quote=F,row.names=F,col.names=F,sep=' ',append=T)

write.table('#',paste(nam_name,'_reseq_weighted_genos_trim_yeshet_assigned_genos.txt',sep=''),
            quote=F,row.names=F,col.names=F,sep=' ',append=T)

print(c(ncol(genos_only),nrow(genos_only),0))

for(i in 1:nrow(out)){
  print(i)
  write.table(t(c(loci[i],paste(out[i,],sep=''))),paste(nam_name,'_reseq_weighted_genos_trim_yeshet_assigned_genos.txt',sep=''),
              quote=F,row.names=F,col.names=F,sep='',append=T)  
}
out2=cbind(genos[,1:2],assigned)
colnames(out2)=c(colnames(genos[,1:2]),colnames(genos[,10:ncol(genos)]))

write.table(out2,paste(nam_name,'_reseq_weighted_genos_trim_yeshet_assigned_genos.tsv',sep=''),
            quote=F,row.names=F,col.names=T,sep='\t')  


