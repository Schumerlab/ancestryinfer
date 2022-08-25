arrArgs <- commandArgs(trailingOnly = TRUE);

if(length(arrArgs)<2){
stop("usage is: Rscript recode_bed_files.R bed1 bed2\n") 
}

bed1<-as.character(arrArgs[1])
bed2<-as.character(arrArgs[2])

out<-paste(bed1,"_combined",sep="")
final<-paste(bed1,"_combined_recode",sep="")

bedcommand<-paste("bedtools intersect -a ",bed1," -b ",bed2," -wo | sort -n -k2,2 -k7,7 > ",out,sep="")

system(bedcommand)

data<-read.csv(file=out,sep="\t",head=FALSE)


recode<-{}
start=1; stop=0;
for(x in 1:length(data$V1)){
hap1=data$V4[x]; hap2=data$V9[x];

if((hap1 == 0) & (hap2 == 0)){geno=0}
if((hap1 == 1) & (hap2 == 1)){geno=2}
if((hap1 == 2) & (hap2 == 2)){geno=5}

if((hap1 == 0) & (hap2 == 2)){geno=3}
if((hap1 == 2) & (hap2 == 0)){geno=3}

if((hap1 == 2) & (hap2 == 1)){geno=4}
if((hap1 == 1) & (hap2 == 2)){geno=4}

if((hap1 == 0) & (hap2 == 1)){geno=1}
if((hap1 == 1) & (hap2 == 0)){geno=1}

chrom=data[x,][1]
start=stop+1; stop=data$V11[x]+start
recode<-rbind(recode,cbind(chrom,start,stop,geno))
}

write.table(recode,file=final,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
