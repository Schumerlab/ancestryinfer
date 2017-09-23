#perl! -w

if(@ARGV<2){
    print "perl collapse_parental_counts_1SNPfilter_v2.pl combined_infile read_length\n";
}

my $infile=shift(@ARGV); chomp $infile;
open IN, $infile or die "cannot open par1 mapped infile\n";

my $read_length=shift(@ARGV); chomp $read_length;

open OUT1, ">"."$infile".".pass.formatted";
open OUT2, ">"."$infile"."parental.format";
open FAIL, ">"."$infile".".failingsites";

my $par1count_map1=0; my $par2count_map1=0;
my $par1count_map2=0; my $par2count_map2=0;
my @failed_sites;
my $focal_site="";
my $par1_ref="";
my $par1_alt="";
my $par2_ref="";
my $par2_alt="";
my $recrate="";
my $pos_prev=1;
my $distance=0;

while(my $line=<IN>){

    chomp $line;

    my @data=split(/\t/,$line);

    $focal_site=$data[0]; chomp $focal_site;
    (my $group, my $site)=split(/\_/,$focal_site);
    chomp $site; chomp $group;

    $par1_ref=$data[1]; chomp $par1_ref;
    $par1_alt=$data[2]; chomp $par1_alt;

    $par2_ref=$data[3]; chomp $par2_ref;
    $par2_alt=$data[4]; chomp $par2_alt;

    $recrate=$data[5]; chomp $recrate;

    $par1count_map1=$data[6]; chomp $par1count_map1;
    $par2count_map1=$data[7]; chomp $par2count_map1;

    $par1count_map2=$data[15]; chomp $par1count_map2;
    $par2count_map2=$data[14]; chomp $par2count_map2;

    my $sum1=$par1count_map1+$par2count_map1;
    my $sum2=$par1count_map2+$par1count_map2;

    $distance=$site-$pos_prev;
    $pos_prev=$site;

    if($distance >= $read_length){

    if(($sum1 != 0) && ($sum2 == 0)){
       #print "$sum1\t$sum2\n";
	print OUT2 "$group\t$site\t$par1_ref\t$par1_alt\t$par2_ref\t$par2_alt\t$recrate\n";
	print OUT1 "$par1count_map1\t$par2count_map1\n";
       }#mapping to parent2 gives no info
    if(($sum1 == 0) && ($sum2 != 0)){
	print OUT2 "$group\t$site\t$par1_ref\t$par1_alt\t$par2_ref\t$par2_alt\t$recrate\n";
	print OUT1 "$par1count_map2\t$par2count_map2\n";
      }#mapping to parent1 gives no info
    if(($sum1 != 0) && ($sum2 != 0)){
	
	if(($par1count_map1 == $par1count_map2) && ($par2count_map1 == $par2count_map2)){

	    print OUT2 "$group\t$site\t$par1_ref\t$par1_alt\t$par2_ref\t$par2_alt\t$recrate\n";
	    print OUT1 "$par1count_map1\t$par2count_map1\n";

	}#if they agree this is a passing base
	else{

            print OUT2 "$group\t$site\t$par1_ref\t$par1_alt\t$par2_ref\t$par2_alt\t$recrate\n";
            print OUT1 "0\t0\n";
            print FAIL "FAILING\t$par1count_map1\t$par2count_map1\t$par1count_map2\t$par2count_map2\n";

	}#this is a failing site, catalog

    }#both non-zero, check to see if base is passing (in agreement)
    
    if(($sum1 == 0) && ($sum2 == 0)){

	    print OUT2 "$group\t$site\t$par1_ref\t$par1_alt\t$par2_ref\t$par2_alt\t$recrate\n";
	    print OUT1 "0\t0\n";
	    print FAIL "FAILING\t$par1count_map1\t$par2count_map1\t$par1count_map2\t$par2count_map2\n";

	}#this is a failing site, catalog 

    }#passes 1 SNP per read threshold
    else{

	print OUT2 "$group\t$site\t$par1_ref\t$par1_alt\t$par2_ref\t$par2_alt\t$recrate\n";
	print OUT1 "0\t0\n";
	print FAIL "MASKED\t$par1count_map1\t$par2count_map1\t$par1count_map2\t$par2count_map2\n";

    }#treat as missing site   

}#for all lines

