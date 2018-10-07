#perl! -w

my $list=shift(@ARGV); chomp $list;
open IN, $list or die "cannot open list of results\n";

my $original_names=shift(@ARGV); chomp $original_names;
open NAMES, $original_names or die "cannot open list of file names\n";

my $save_files=shift(@ARGV); chomp $save_files;

my $tag="";

my $focalchr=shift(@ARGV); chomp $focalchr;
if($focalchr ne 0){
open TAGCHR, $focalchr or die "cannot open $focalchr file\n"
}# if there are focal chrs 
else{ $tag="$tag"."_"."allchrs";}
while(my $chrnames=<TAGCHR>){
    chomp $chrnames;
    $tag="$tag"."_"."$chrnames";
}#add name of focal chroms
#else{ $tag="$tag"."_"."allchrs";}
print "file tag is $tag\n";

my @filenames=();
while(my $n=<NAMES>){

    chomp $n;
    my @subnames=split(/\t/,$n);
    my $nm=$subnames[0]; chomp $nm;
    my @idparsed=split(/\//,$nm);
#    print "@idparsed\n";
    my $length=scalar(@idparsed) -1;
#    print "$length\n";
    my $trimmedname=$idparsed[$length]; chomp $trimmedname;
#    print "$trimmedname\n";
    $trimmedname=~ s/\.gz//g;
    push(@filenames, $trimmedname);

}#collect original file names

my $par1_string="";
my $par2_string="";
my $counter=0;
while(my $line=<IN>){

    $counter=$counter+1;
    chomp $line;
    my @elements=split(/\t/,$line);
    my $focal_file=$elements[0]; chomp $focal_file;
    $focal_file="$focal_file".".posterior";

    #match filename to simpler name
    my $currname="";
    for my $k (0..scalar(@filenames)-1){ 
	my $focalmatch=$filenames[$k]; chomp $focalmatch;

	if($focal_file =~ /$focalmatch/g){
	    print "matching names: $focalmatch\t$focal_file\n";
	    $currname=$focalmatch;
	}

    }#for all filenames

    if($counter==1){
	
    my $temp_par1="$focal_file".".par1.temp";
    system("echo '\t$currname' > $temp_par1");
    system("cut -f 1-2,3 $focal_file | tail -n +2 >> $temp_par1");

    my $temp_par2="$focal_file".".par2.temp";
    system("echo '\t$currname' > $temp_par2");
    system("cut -f 1-2,5 $focal_file | tail -n +2 >> $temp_par2");

    $par1_string="$temp_par1";
    $par2_string="$temp_par2";

    system("perl -pi -e 's/^(.+?)\t/\\1:/g' $temp_par1");
    system("perl -pi -e 's/^(.+?)\t/\\1:/g' $temp_par2");

    }#generate header information
    else{
    my $temp_par1="$focal_file".".par1.temp";
    system("echo '$currname' > $temp_par1");
    system("cut -f 3 $focal_file | tail -n +2 >> $temp_par1"); 
    
    my $temp_par2="$focal_file".".par2.temp";
    system("echo '$currname' > $temp_par2");
    system("cut -f 5 $focal_file | tail -n +2 >> $temp_par2");

    $par1_string="$par1_string"." "."$temp_par1";
    $par2_string="$par2_string"." "."$temp_par2";

    }

    if($save_files==0){
    system("rm $focal_file"); 
    }#remove individual posterior probs files

}#for all individuals run through the HMM

my $file1name="ancestry-probs-par1_transposed"."$tag".".tsv";
my $file2name="ancestry-probs-par2_transposed"."$tag".".tsv";
system("paste $par1_string > $file1name");
system("paste $par2_string > $file2name");
if($save_files==0){
system("rm $par1_string");
system("rm $par2_string");
system("rm map_batch*.sh samtools_batch*.sh split_jobs_list sam_files_mapped_to_parent1 sam_files_mapped_to_parent2 hmm_batch.sh");
#print "$par1_string\n";

system("rm HMM.parental.files.list*");
system("rm HMM.hybrid.files.list*");
#my $prefix="$original_names".".";
#system("rm $prefix*");
system("rm current.samples.list");

}#remove files, don't save
