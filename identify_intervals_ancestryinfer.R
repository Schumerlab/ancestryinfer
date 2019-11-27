arrArgs<-commandArgs(trailingOnly = TRUE);
genos<-as.character(arrArgs[1])
path<-as.character(arrArgs[2])

if(length(arrArgs)<2){
stop("usage is: Rscript identify_intervals_ancestryinfer.R genotypes_file_name path_to:transpose_nameout.pl\n");
}#print usage

command1<-paste("perl ",path,"/transpose_nameout.pl ",genos,sep="")
system(command1)

genos_transposed<-paste(genos,"_transposed",sep="")

select_markers<-paste("cut -f 1 ",genos_transposed," | perl -p -e 's/:/\t/g' | cut -f 1 | uniq | tail -n +2 > ",genos,"_header",sep="")
system(select_markers)

chroms_raw<-read.csv(file=paste(genos,"_header",sep=""),sep="\t",head=FALSE)
chroms<-chroms_raw$V1

outfile<-paste(genos,"_ancestrytransitions_allchrs",sep="")
file.remove(outfile)

for(k in 1:length(chroms)){

command2=paste("grep ",chroms[k]," ",genos_transposed," | perl -p -e ","'","s/:/\t/g","'"," > ",genos_transposed,"_",chroms[k],sep="")
system(command2)

file=paste(genos_transposed,"_",chroms[k],sep="")
data<-read.csv(file=file,sep="\t",head=FALSE,as.is=T)

intervals<-{}
count=0

for(y in 3:(length(data[1,]))){

focal_geno=data[,y]
sites<-data[,2]
na_count=0;
geno_prev=subset(focal_geno,!is.na(focal_geno)==TRUE)[1];
first_region=1
last_geno=focal_geno[1]
start=sites[1]

for(x in 1:length(focal_geno)){
geno_current=focal_geno[x];

if(!is.na(geno_current)==TRUE){

	if(geno_current == geno_prev){

	} else{

	if(is.na(last_geno)==TRUE){
	stop=sites[x]
	geno_prev=geno_current

	intervals<-rbind(intervals,cbind(start,stop,y-2))

	}#make sure previous genotype was NA for this mode

	}#geno_current does not equal geno_previous

} else{

if((is.na(geno_current)==TRUE)	&(is.na(last_geno)==FALSE)){
	start=sites[x-1]
}#option 1: this is the first site in the interval region

}#if is/is not NA


if((is.na(geno_current)==FALSE)&(is.na(last_geno)==FALSE)&(geno_current != last_geno)){
	start=sites[x-1]
	stop=sites[x]
	geno_prev=geno_current

	intervals<-rbind(intervals,cbind(start,stop,y-2))
	count=count+1
}#option 2: the interval was between this marker and the last marker


last_geno=geno_current

}#for all sites in the focal individual

}#for all individuals

write.table(cbind(rep(as.character(chroms[k]),length(intervals[,1])),intervals),file=outfile,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
#write.table(intervals,file=outfile,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)

command3<-paste(genos,"_",chroms[k],sep="")
#system(command3)

}#all chroms

