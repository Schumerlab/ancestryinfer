#perl! -w

my $list=shift(@ARGV); chomp $list;

my $original_names=shift(@ARGV); chomp $original_names;

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

$list="$list"."$tag";
open IN, $list or die "cannot open $list list of results\n";

$original_names="$original_names"."$tag";
open NAMES, $original_names or die "cannot open $original_names list of file names\n";

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
my $par3_string="";
my $het1_string="";
my $het2_string="";
my $het3_string="";
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
	$focalmatch =~ s/$tag//g;
        #print "$focalmatch\n";
	if($focal_file =~ /$focalmatch/g){
	    print "matching names: $focalmatch\t$focal_file\n";
	    $currname=$focalmatch;
	}

    }#for all filenames

    if($counter==1){
	
    my $temp_par1="$focal_file".".par1.results";
    system("echo '\t$currname' > $temp_par1");
    system("cut -f 1-2,3 $focal_file | tail -n +2 >> $temp_par1");

    my $temp_par2="$focal_file".".par2.results";
    system("echo '\t$currname' > $temp_par2");
    system("cut -f 1-2,6 $focal_file | tail -n +2 >> $temp_par2");

    my $temp_par3="$focal_file".".par3.results";
    system("echo '\t$currname' > $temp_par3");
    system("cut -f 1-2,8 $focal_file | tail -n +2 >> $temp_par3");

    my $temp_het1="$focal_file".".par1par2.results";
    system("echo '\t$currname' > $temp_het1");
    system("cut -f 1-2,4 $focal_file | tail -n +2 >> $temp_het1");

    my $temp_het2="$focal_file".".par1par3.results";
    system("echo '\t$currname' > $temp_het2");
    system("cut -f 1-2,5 $focal_file | tail -n +2 >> $temp_het2");

    my $temp_het3="$focal_file".".par2par3.results";
    system("echo '\t$currname' > $temp_het3");
    system("cut -f 1-2,7 $focal_file | tail -n +2 >> $temp_het3");

    $par1_string="$temp_par1";
    $par2_string="$temp_par2";
    $par3_string="$temp_par3";
    $het1_string="$temp_het1";
    $het2_string="$temp_het2";
    $het3_string="$temp_het3";

    system("perl -pi -e 's/^(.+?)\t/\\1:/g' $temp_par1");
    system("perl -pi -e 's/^(.+?)\t/\\1:/g' $temp_par2");
    system("perl -pi -e 's/^(.+?)\t/\\1:/g' $temp_par3");
    system("perl -pi -e 's/^(.+?)\t/\\1:/g' $temp_het1");
    system("perl -pi -e 's/^(.+?)\t/\\1:/g' $temp_het2");
    system("perl -pi -e 's/^(.+?)\t/\\1:/g' $temp_het3");

    }#generate header information
    else{
    my $temp_par1="$focal_file".".par1.results";
    system("echo '$currname' > $temp_par1");
    system("cut -f 3 $focal_file | tail -n +2 >> $temp_par1"); 
    
    my $temp_par2="$focal_file".".par2.results";
    system("echo '$currname' > $temp_par2");
    system("cut -f 6 $focal_file | tail -n +2 >> $temp_par2");

    my $temp_par3="$focal_file".".par3.results";
    system("echo '\t$currname' > $temp_par3");
    system("cut -f 8 $focal_file | tail -n +2 >> $temp_par3");

    my $temp_het1="$focal_file".".par1par2.results";
    system("echo '\t$currname' > $temp_het1");
    system("cut -f 4 $focal_file | tail -n +2 >> $temp_het1");

    my $temp_het2="$focal_file".".par1par3.results";
    system("echo '\t$currname' > $temp_het2");
    system("cut -f 5 $focal_file | tail -n +2 >> $temp_het2");

    my $temp_het3="$focal_file".".par2par3.results";
    system("echo '\t$currname' > $temp_het3");
    system("cut -f 7 $focal_file | tail -n +2 >> $temp_het3");

    $par1_string="$par1_string"." "."$temp_par1";
    $par2_string="$par2_string"." "."$temp_par2";
    $par3_string="$par3_string"." "."$temp_par3";
    $het1_string="$het1_string"." "."$temp_het1";
    $het2_string="$het2_string"." "."$temp_het2";
    $het3_string="$het3_string"." "."$temp_het3";
    }

    if($save_files==0){
    system("rm $focal_file"); 
    }#remove individual posterior probs files

}#for all individuals run through the HMM

my $file1name="ancestry-probs-par1_transposed"."$tag".".tsv";
my $file2name="ancestry-probs-par2_transposed"."$tag".".tsv";
my $file3name="ancestry-probs-par3_transposed"."$tag".".tsv";
my $file4name="ancestry-probs-par1par2_transposed"."$tag".".tsv";
my $file5name="ancestry-probs-par1par3_transposed"."$tag".".tsv";
my $file6name="ancestry-probs-par2par3_transposed"."$tag".".tsv";
system("paste $par1_string > $file1name");
system("paste $par2_string > $file2name");
system("paste $par3_string > $file3name");
system("paste $het1_string > $file4name");
system("paste $het2_string > $file5name");
system("paste $het3_string > $file6name");

if($save_files==0){
system("rm $par1_string");
system("rm $par2_string");
system("rm $par3_string");
system("rm $het1_string");
system("rm $het2_string");
system("rm $het3_string");
system("rm map_batch*.sh samtools_batch*.sh sam_files_mapped_to_parent1 sam_files_mapped_to_parent2");
#print "$par1_string\n";

}#remove files, don't save
