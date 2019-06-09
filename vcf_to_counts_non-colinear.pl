#perl! -w

my $infile=shift(@ARGV); chomp $infile;

my $aims_list=shift(@ARGV); chomp $aims_list;

my $path=shift(@ARGV); chomp $path;

my $outfile="$infile"."_counts";
open OUT, ">$outfile";

####generate bed file and overlap with input vcf
my $aims_bed="$aims_list".".mod";
my $infile_aims="$infile".".aims";

my $vcf_mod="$infile".".mod";
system("cat $infile | perl -p -e 's/_/\t/g' | awk -v OFS=\'\\t\' \'\$1=\$1\"\_\"\$2\' > $vcf_mod");
system("perl $path/combine_FAS_scriptome_snippet.pl $aims_bed $vcf_mod $infile_aims");
open IN, $infile_aims or die "cannot open modified mpileup"; 

my $compound=""; my $group=""; my $pos=""; my $reads1=0; my $reads2=0; my $ref=""; my $alt=""; my $match="";
while(my $line = <IN>){

    chomp $line;
    
    if(($line !~ /\#/g) && ($line !~ /INDEL/g)){
	
	$line =~ s/_/\t/g;
	#print "$line\n"; #check filter
	my @elements=split(/\t/,$line);
	$group=$elements[0]; chomp $group;
	$pos=$elements[1]; chomp $pos;
	$ref=$elements[3]; chomp $ref;
	$alt=$elements[4]; chomp $alt;
	#print "$group\t$pos\t$ref\t$alt\n"; #check assignment

	$match=$elements[7]; chomp $match;

	if($match eq $pos){

#####change here is any change in mpileup format occurs	
	my @info=split(/\=/,$elements[13]);
	#print "@info\n";

	my $genostring=""; my @genodat=();
	
	$genostring=$elements[15]; chomp $genostring;
	@genodat=split(/:/,$genostring);

	my $geno=$genodat[0];
	$geno =~ s/\//_/g;
	#print "$geno\n";

	my $col=scalar(@info)-2;
	#print "$col\n";

	if($geno eq '1_1'){
	my @remainder=split(/;/,$info[$col]);

	#print "$remainder[0]","\t","$info[2]","\n";
	my @deptharray=split(/\,/,$remainder[0]); 
	#print "@deptharray\n";
	my $refcount=$deptharray[0]+$deptharray[1];
	my $altcount=$deptharray[2]+$deptharray[3];

	my @bparray=();
	push(@bparray,($ref,$alt));
	my @sortedbp=sort { lc($a) cmp lc($b) } @bparray;

	#check that ref and alt are the same
	my $ref_alt="$elements[9]"."$elements[10]"; my $alt_ref="$elements[10]"."$elements[9]"; my $true_ref_alt="$ref"."$alt";
	#print "$elements[9]\t$elements[10]\n";
	if(($ref_alt eq $true_ref_alt) or ($alt_ref eq $true_ref_alt)){
	print OUT "$group"."_"."$pos"."_"."$sortedbp[0]"."_"."$sortedbp[1]\t$refcount\t$altcount\n";
	}#check alleleic concordance

	#print "$geno\t$line\n";
	}#homo
	elsif($geno eq '0_1'){

	    my @remainder=split(/;/,$info[$col]);
	    #print "$remainder[0]","\t","$info[2]","\n";                                                      
	    my @deptharray=split(/\,/,$remainder[0]);
	    #print "@deptharray\n";                                                                                                                                                   
	    my $refcount=$deptharray[0]+$deptharray[1];
	    my $altcount=$deptharray[2]+$deptharray[3];

	    my @bparray=();
	    push(@bparray,($ref,$alt));
	    my @sortedbp=sort { lc($a) cmp lc($b) } @bparray;
	    #check that ref and alt are the same                                                                             
            my $ref_alt="$elements[9]"."$elements[10]"; my $alt_ref="$elements[10]"."$elements[9]"; my $true_ref_alt="$ref"."$alt";
	    if(($ref_alt eq$true_ref_alt) or ($alt_ref eq $true_ref_alt)){
	    print OUT "$group"."_"."$pos"."_"."$sortedbp[0]"."_"."$sortedbp[1]\t$refcount\t$altcount\n";
	    }#check alleleic concordance    
	}#hets
	elsif($geno eq '0_0'){

	    my @remainder=split(/;/,$info[$col]);
            #print "$remainder[0]","\t","$info[2]","\n";
            my @deptharray=split(/\,/,$remainder[0]);
            #print "@deptharray\n";
            my $refcount=$deptharray[0]+$deptharray[1];
            my $altcount=$deptharray[2]+$deptharray[3];

	    my @bparray=();
            push(@bparray,($ref,$alt));
            my @sortedbp=sort { lc($a) cmp lc($b) } @bparray;
            print OUT "$group"."_"."$pos"."_"."$sortedbp[0]"."_"."$sortedbp[1]\t$refcount\t$altcount\n";

	}
	#print "$refcount\t$altcount\n";

	}#Ensure perfect overlap from bedtools

    }#if not a header line


}#all lines
