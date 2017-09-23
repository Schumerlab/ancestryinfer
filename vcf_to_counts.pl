#perl! -w

my $infile=shift(@ARGV); chomp $infile;
open IN, $infile or die "cannot open mpileup\n";

my $outfile="$infile"."_counts";
open OUT, ">$outfile";

my $group=""; my $pos=""; my $reads1=0; my $reads2=0; my $ref=""; my $alt="";
while(my $line = <IN>){

    chomp $line;

    if(($line !~ /\#/g) && ($line !~ /INDEL/g)){

	#print "$line\n"; #check filter
	my @elements=split(/\t/,$line);
	$group=$elements[0]; chomp $group;
	$pos=$elements[1]; chomp $pos;
	$ref=$elements[3]; chomp $ref;

	$alt=$elements[4]; chomp $alt;
	
	my @info=split(/\=/,$elements[7]); #change here if change mpileup format in any way
	
	my $genostring=$elements[9]; chomp $genostring;
	my @genodat=split(/:/,$genostring);
	
	my $geno=$genodat[0];
	$geno =~ s/\//_/g;

	my $col=scalar(@info)-2;
	#print "$col\n";

	if($geno eq '1_1'){
	my @remainder=split(/;/,$info[$col]);

	#print "$remainder[0]","\t","$info[2]","\n";
	my @deptharray=split(/\,/,$remainder[0]); 
	#print "@deptharray\n";
	my $refcount=$deptharray[0]+$deptharray[1];
	my $altcount=$deptharray[2]+$deptharray[3];

	print OUT "$group"."_"."$pos\t$ref\t$alt\t$refcount\t$altcount\n";

	#print "$geno\t$line\n";
	}#homo
	elsif($geno eq '0_1'){

	    my @remainder=split(/;/,$info[$col]);

	    #print "$remainder[0]","\t","$info[2]","\n";                                                                                                                              
	    my @deptharray=split(/\,/,$remainder[0]);
        #print "@deptharray\n";                                                                                                                                                   
	    my $refcount=$deptharray[0]+$deptharray[1];
	    my $altcount=$deptharray[2]+$deptharray[3];

	    print OUT "$group"."_"."$pos\t$ref\t$alt\t$refcount\t$altcount\n";

	}#hets
	else{

	    my @doublealt=split(/,/,$alt);
	    my $alt1=$doublealt[0]; chomp $alt1;
	    my $alt2=$doublealt[1]; chomp $alt2;
	    my @remainder=split(/;/,$info[$col]);

            #print "$remainder[0]","\t","$info[2]","\n";
            my @deptharray=split(/\,/,$remainder[0]);
            #print "@deptharray\n";
            my $refcount=$deptharray[0]+$deptharray[1];
            my $altcount=$deptharray[2]+$deptharray[3];

            print OUT "$group"."_"."$pos\t$alt1\t$alt2\t$refcount\t$altcount\n";

	}
	#print "$refcount\t$altcount\n";

    }#if not a header line


}#all lines
