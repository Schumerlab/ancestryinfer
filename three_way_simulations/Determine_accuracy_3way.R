arrArgs <- commandArgs(trailingOnly = TRUE);

if(length(arrArgs)<4){
stop("usage is: Rscript Determine_accuracy.R true_ancestry.bed individ genotypes_file path\n") 
}

binsfile<-as.character(arrArgs[1])
indiv<-as.character(arrArgs[2])
infile<-as.character(arrArgs[3])
#pp1<-as.character(arrArgs[4])
#pp2<-as.character(arrArgs[5])
path<-as.character(arrArgs[4])

bins<-read.csv(file=binsfile,sep="\t",head=FALSE)


bins<-bins[order(bins$V2),]
#command 0 set
#head(bins)

command0=paste("head -n 1 ",infile," > accuracy_",indiv,"_","genotypes_file",sep="")

system(command0)

#command0b=paste("head -n 1 ",infile," > accuracy_",indiv,"_","genotypes_file_pp1",sep="")
#command0c=paste("head -n 1 ",infile," > accuracy_",indiv,"_","genotypes_file_pp2",sep="")

#system(command0b)
#system(command0c)

#command 1 set
command1=paste("grep ",indiv,"_ ",infile," ",">> accuracy_",indiv,"_","genotypes_file",sep="")

system(command1)

#command1b=paste("grep ",indiv,"_ ",pp1," ",">> accuracy_",indiv,"_","genotypes_file_pp1",sep="")
#command1c=paste("grep ",indiv,"_ ",pp2," ",">> accuracy_",indiv,"_","genotypes_file_pp2",sep="")

#system(command1b)
#system(command1c)

command2=paste("perl ",path,"/transpose_genotypes_tsv.pl ","accuracy_",indiv,"_","genotypes_file",sep="")
system(command2)

#command2b=paste("perl ",path,"/transpose_genotypes_tsv.pl ","accuracy_",indiv,"_","genotypes_file_pp1",sep="")
#command2c=paste("perl ",path,"/transpose_genotypes_tsv.pl ","accuracy_",indiv,"_","genotypes_file_pp2",sep="")

#system(command2b)
#system(command2c)

#command 2 set
data<-read.csv(file=paste("accuracy_",indiv,"_","genotypes_file_transposed",sep=""),sep="\t",head=TRUE)
#pp1<-read.csv(file=paste("accuracy_",indiv,"_","genotypes_file_pp1_transposed",sep=""),sep="\t",head=TRUE)
#pp2<-read.csv(file=paste("accuracy_",indiv,"_","genotypes_file_pp2_transposed",sep=""),sep="\t",head=TRUE)

options(scipen=999)
whole_genome<-{}
start=1
stop=bins$V3[1]
counts_het1=0
counts_het2=0
counts_het3=0
counts_par1=0
counts_par2=0
counts_par3=0
accurate_counts=0
inaccurate_counts=0
mean_pos<-{}
for (x in 1:length(bins$V1)){

focal<-subset(data,data$pos>=start & data$pos<=stop)

counts_par1=length(subset(focal[,1],focal[,3] ==0)) #par1
counts_par2=length(subset(focal[,1],focal[,3] ==2)) #par2
counts_par3=length(subset(focal[,1],focal[,3] ==5)) #par3

counts_het1=length(subset(focal[,1],focal[,3] ==1)) #par1par2
counts_het2=length(subset(focal[,1],focal[,3] ==3)) #par1par3
counts_het3=length(subset(focal[,1],focal[,3] ==4)) #par2par3

if(bins$V4[x] == 1){
accurate_counts=accurate_counts+counts_het1
inaccurate_counts=inaccurate_counts+0.5*counts_par1+0.5*counts_par2+counts_par3+0.5*counts_het2+0.5*counts_het3
}
if(bins$V4[x] == 3){
accurate_counts=accurate_counts+counts_het2
inaccurate_counts=inaccurate_counts+0.5*counts_par1+counts_par2+0.5*counts_par3+0.5*counts_het1+0.5*counts_het3
}
if(bins$V4[x] == 4){
accurate_counts=accurate_counts+counts_het3
inaccurate_counts=inaccurate_counts+counts_par1+0.5*counts_par2+0.5*counts_par3+0.5*counts_het1+0.5*counts_het2
}


if(bins$V4[x] == 0){
accurate_counts=accurate_counts+counts_par1
inaccurate_counts=inaccurate_counts+0.5*counts_het1+0.5*counts_het2+counts_het3+counts_par2+counts_par3
}
if(bins$V4[x] == 2){
accurate_counts=accurate_counts+counts_par2
inaccurate_counts=inaccurate_counts+0.5*counts_het1+counts_het2+0.5*counts_het3+counts_par1+counts_par3
}
if(bins$V4[x] == 5){
accurate_counts=accurate_counts+counts_par3
inaccurate_counts=inaccurate_counts+counts_het1+0.5*counts_het2+0.5*counts_het3+counts_par1+counts_par2
}


#write.table(cbind(mean_pos))

#head(cbind(indiv,start,stop,counts_het,counts_par1,counts_par2,bins$V4[x],accurate_counts,inaccurate_counts,mean(mean_pos)))

whole_genome<-rbind(whole_genome,cbind(indiv,start,stop,bins$V4[x],counts_het1,counts_het2,counts_het3,counts_par1,counts_par2,counts_par3,accurate_counts,inaccurate_counts,accurate_counts/(accurate_counts+inaccurate_counts)))

counts_het1=0
counts_het2=0
counts_het3=0
counts_par1=0
counts_par2=0
counts_par3=0

start=bins$V2[x+1]
stop=bins$V3[x+1]
accurate_counts=0
inaccurate_counts=0
}

write.table(whole_genome,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
