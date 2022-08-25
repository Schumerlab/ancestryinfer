#perl! -w

if(@ARGV<9){
    print "perl simulate_3way_hybrids.pl genome1 genome2 genome3 SELAM_tracts_file chrom_length_bp rec_rate_M_bp focal_chr num_reads simulation_name\n"; exit;
}

my $genome1=shift(@ARGV); chomp $genome1;
my $genome2=shift(@ARGV); chomp $genome2;
my $genome3=shift(@ARGV); chomp $genome3;

my $tracts=shift(@ARGV); chomp $tracts;
open IN, $tracts or die "cannot open SELAM tracts file\n";

my $chromlength=shift(@ARGV); chomp $chromlength;

my $rec=shift(@ARGV); chomp $rec;

my $chr=shift(@ARGV); chomp $chr;

my $num_reads=shift(@ARGV); chomp $num_reads;

my $out=shift(@ARGV); chomp $out;

my $bed_list="$out"."_bed_list";
open BED, ">$bed_list";

my $curr_indiv=""; my $curr_hap=0; my $prev_indiv=""; 
my $file=""; my $tracts=""; my $hap1=""; my $hap2=""; my $current1=0; my $current2=0; my $total=$chromlength;
my $start1=0; my $stop1=0;
my $start2=0; my $stop2=0;
my $hap_file=""; my $tracts_file="";
my $add1=""; my $add2="";
my $sex="";

my $counter=0;
while (my $line=<IN>){

    $counter++;

    if($line !~ /#/g){
    chomp $line;
    my @elements=split(/\t/,$line);
    my $curr_indiv_tmp=$elements[3];
    $sex=$elements[2];
    $curr_indiv="$curr_indiv_tmp"."_sex"."$sex";
    $curr_hap=$elements[5];
    $draw1=$elements[6];
    $hap_start=$elements[7]; $hap_stop=$elements[8];
    my $hap_length=$hap_stop-$hap_start;
    
    my $tract=int($hap_length*(1/$rec));

#    print "tract length is $tract\n";
    if($counter eq 1){
        $file="$out"."_$curr_indiv".".fa";
        $tracts1="$out"."_$curr_indiv"."_tracts_hap1.bed";
        $tracts2="$out"."_$curr_indiv"."_tracts_hap2.bed";
        open OUT, ">$file";
        open TRACTS1, ">$tracts1";
        open TRACTS2, ">$tracts2";
	print BED "$tracts1\t$tracts2\n";
	$prev_indiv=$curr_indiv;
    }#first individual   



    if($curr_indiv ne $prev_indiv){

	print "previous individual was $prev_indiv, current is $curr_indiv\n";

        print OUT ">"."indiv"."$prev_indiv"."_haps\n"."$hap1"."$hap2\n";
	#print OUT ">"."indiv"."$prev_indiv"."_hap2\n"."$hap2\n";

	my $current="$out"."_$prev_indiv".".fa";
	my $read1="$out"."_$prev_indiv"."_r1.fq";
	my $read2="$out"."_$prev_indiv"."_r2.fq";
	my $log="log_indiv"."$prev_indiv";

	system("wgsim -N $num_reads -1 150 -2 150 -S 1 -e 0.01 -r 0.05 -R 1 $current $read1 $read2 > $log");
	print "zipping reads $read1 $read2\n";
	system("gzip $read1 $read2");

        $file="$out"."_$curr_indiv".".fa";
	$tracts1="$out"."_$curr_indiv"."_tracts_hap1.bed";
	$tracts2="$out"."_$curr_indiv"."_tracts_hap2.bed";
	open OUT, ">$file";
	open TRACTS1, ">$tracts1";
	open TRACTS2, ">$tracts2";
	print BED "$tracts1\t$tracts2\n";

    $hap1=""; $hap2="";
    $current1=0; $current2=0; $total=$chromlength;
    $start1=0; $stop1=0;
    $start2=0; $stop2=0;
    $add1=""; $add2="";
}  #clear between individuals

    $prev_indiv=$curr_indiv;

	if ($draw1 eq 0){ $hap_file = $genome1;}
	elsif ($draw1 eq 1){ $hap_file = $genome2;}
	elsif ($draw1 eq 2){ $hap_file = $genome3;}

    if($curr_hap eq 0){
	$start1=$stop1+1; $stop1=$start1+$tract-1;
        $current1=$current1+$tract;
        
	#print "hap1\t$hap_length\t$start1\t$stop1\t$curr_indiv\t$prev_indiv\t$current1\n";

	if ($current1>=$chromlength){$stop1=$chromlength; $current1=$chromlength;}
	if($start1 < $chromlength){
	$add1=qx(fastahack $hap_file -r $chr:$start1..$stop1); chomp $add1;
	print TRACTS1 "$chr\t$start1\t$stop1\t$draw1\tindiv"."$curr_indiv"."_hap1\n";
	$hap1="$hap1"."$add1";
	}#fill only to the chromosome edge then stop
	
	my $test_length=length($hap1);
	print "hap1\t$hap_length\t$start1\t$stop1\t$curr_indiv\t$prev_indiv\t$current1\t$test_length\n";

    }
    if($curr_hap eq 1){
	$start2=$stop2+1; $stop2=$start2+$tract-1;
        $current2=$current2+$tract;

	#print "hap2\t$hap_length\t$start2\t$stop2\t$curr_indiv\t$prev_indiv\t$current2\n";

        if ($current2>=$chromlength){$stop2=$chromlength; $current2=$chromlength;}
	if($start2 < $chromlength){
        $add2=qx(fastahack $hap_file -r $chr:$start2..$stop2); chomp $add2;
	print TRACTS2 "$chr\t$start2\t$stop2\t$draw1\tindiv"."$curr_indiv"."_hap2\n";
	$hap2="$hap2"."$add2";
	}#fill only to the chromosome edge then stop 
	
	my $test_length=length($hap2);
	print "hap2\t$hap_length\t$start2\t$stop2\t$curr_indiv\t$prev_indiv\t$current2\t$test_length\n";

    }

    }#for all data lines

}#generate data for all individuals


print OUT ">"."indiv"."$curr_indiv"."_haps\n"."$hap1"."$hap2\n";
#print OUT ">"."indiv"."$curr_indiv"."_hap2\n"."$hap2\n";

my $current="$out"."_$curr_indiv".".fa";
my $read1="$out"."_$curr_indiv"."_r1.fq";
my $read2="$out"."_$curr_indiv"."_r2.fq";
my $log="log_indiv"."$curr_indiv";

system("wgsim -N $num_reads -1 150 -2 150 -S 1 -e 0.01 -r 0.05 -R 1 $current $read1 $read2 > $log");
print "zipping reads\n";
system("gzip $read1 $read2");
