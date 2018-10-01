#perl! -w

if(@ARGV<2){
    print "perl combine_all_individuals_hmm.pl parental_list individual_list parental_prior counts initial_time_admix focal_chrom read_length\n";
}

my $parental=shift(@ARGV); chomp $parental;
open PAR, $parental or die "cannot open parental list file\n";

my $list=shift(@ARGV); chomp $list;
open LIST, $list or die "cannot open individual list file\n";

my $par1_prior=shift(@ARGV); chomp $par1_prior;
if(length($par1_prior)==0){ die "parental priors required\n";}

my $provide_counts=shift(@ARGV); chomp $provide_counts;

my $initial_admix=shift(@ARGV); chomp $initial_admix;

my $focal_chrom=shift(@ARGV); chomp $focal_chrom;

my $read_length=shift(@ARGV); chomp $read_length;

open READLIST, ">current.samples.read.list";

my $prior_admix=0;
if($initial_admix > 0){
    print "prior admixture time provided: $initial_admix\n";
    $prior_admix=1;
}#check prior status

my $set_focal_chrom=0;
if($focal_chrom ne '0'){
    $set_focal_chrom=1;
    print "focal chrom is $focal_chrom\n";
}#focal chrom list provided

    my $aims_list="current_aims_file"; chomp $aims_list;
    open AIMSLIST, $aims_list or die "cannot open AIMs list\n";

    my $aims=<AIMSLIST>; chomp $aims;

    my $aims_length=qx(wc -l $aims | perl -pi -e 's/ +/\t/g' | cut -f 1); chomp $aims_length;

if(-e $provide_counts){

    print "parental counts were provided: $provide_counts\n";
    my $count_length=qx(wc -l $provide_counts | perl -pi -e 's/ +/\t/g' | cut -f 1); chomp $count_length;

    if($aims_length != $count_length){
	die "aims list and parental counts lists have different length, cannot proceed\n";
    }#check

}#check that provided counts are the same length as the aims

open OUT, ">current.samples.list";

my $counter=0;
my $string="";
while ((my $line1=<LIST>) && (my $line2=<PAR>)){

    $counter++;
    chomp $line1; chomp $line2;
    
    my $current_counts=$line1; $current_counts=~ s/\.sam.hmm.combined.pass.formatted//g;
    my $counts_current=qx(wc -l $line1 | perl -p -e 's/ +/\t/g' | cut -f 1); chomp $counts_current;

    if($counts_current eq $aims_length){

    if($counter==1){
	if(length($provide_counts) > 1){
	    $string="$provide_counts"." "."$line1"." ";
	}
	else{
       $string="$line2"." "."$line1"." ";
	}
    } else{
	$string="$string"." "."$line1"." ";
    }

    my $posterior="$line1".".posterior";
    print OUT "$posterior\t2\n";
    print READLIST "$current_counts\n";
    }#ensure proper counts file has been generated, if not warn
    else{print "WARNING $line1 does not contain the appropriate number of markers\n";}

}

print "$string\n";
system("paste $string > all.indivs.hmm.combined");

system("perl -pi -e 's/\\//_/g' current.samples.list");
system("perl -pi -e 's/_/\t/g' all.indivs.hmm.combined");

open CHMM, "all.indivs.hmm.combined" or die "cannot open combined HMM input data\n";
open FILTER, ">all.indivs.hmm.combined.filtered";
my $total=0; my $last_marker=""; my $last_chrom=""; my $curr_marker=""; my $curr_chrom=""; my $rec_distance=0;
####Here, added a way to catalog rec distance for excluded markers
while(my $line=<CHMM>){

    chomp $line;
    my @dataelements=split(/\t/,$line);
    $curr_chrom=$dataelements[0]; chomp $curr_chrom;
    $curr_marker=$dataelements[1]; chomp $curr_marker;
    $total=0;
    for my $i (7..scalar(@dataelements)-1){
        my $focal=$dataelements[$i]; chomp $focal;
        $total=$total+$focal;
    }

    my $distance=$curr_marker-$last_marker; 
    #print "$read_length\t$distance\t$curr_chrom\t$last_chrom\n";
 
   if(($total >0) && (($distance >= $read_length) or ($curr_chrom ne $last_chrom))){
	my $new_line=$line;
    if($rec_distance >0){
	my $replace=$dataelements[6]; chomp $replace;
	#had to modify basic replace fx to something more complicated here due to introduced issue in formatting
	my @splicedline=split(/\t/,$new_line);
	my @trimmedline=@splicedline[ 7 .. $#splicedline ];
	$new_line="$splicedline[0]"."\t"."$splicedline[1]"."\t"."$splicedline[2]"."\t"."$splicedline[3]"."\t"."$splicedline[4]"."\t"."$splicedline[5]"."\t"."$rec_distance"."\t";

	my $end=(scalar(@trimmedline)-1);
	for my $k (0..scalar(@trimmedline)-1){
	    if($k ne $end){
	    $new_line="$new_line"."$trimmedline[$k]"."\t";
	    } else{
	    $new_line="$new_line"."$trimmedline[$k]";
	    }#last line
	}#string out trimmed array
	#print "$new_line\n";
    }#replace recombination distance with current distance
    print FILTER "$new_line\n";
	$rec_distance=0; #reset rec distance tracker
    $last_marker=$curr_marker; $last_chrom=$curr_chrom; #reset marker tracker
    }
    elsif($total==0){
	$rec_distance=$rec_distance+$dataelements[6];
    }#log rec distance of excluded markers
}

####stop here to check for focal chromosome run
my @focal_chrom_array=();
if($set_focal_chrom ==1){
    open FOCAL, "$focal_chrom" or die "cannot open focal chroms file\n";
    while(my $line3=<FOCAL>){
	chomp $line3;
	push(@focal_chrom_array,$line3);
    }#load focal chroms
}#focal chromosomes are defined

###here, subset focal chromosomes to their own file
if($set_focal_chrom ==1){
for my $j (0..scalar(@focal_chrom_array)-1){
    my $focal=$focal_chrom_array[$j];
    print "selecting focal chromosomes $focal\n";
    if($j==0){
	system("grep -w $focal all.indivs.hmm.combined.filtered > all.indivs.hmm.combined.filtered_focalchroms");
    }#first focal chromosome, generate file
    else{
	system("grep -w $focal all.indivs.hmm.combined.filtered >> all.indivs.hmm.combined.filtered_focalchroms");
    }#not first chromosome

  }#grab all focal chromosomes
}#ensure that the switch is flipped before subsetting
if($set_focal_chrom ==0){
system("mv all.indivs.hmm.combined.filtered all.indivs.hmm.combined.filtered_focalchroms"); 
}#all chromosomes are focal, simply rename

####run differently depending on whether prior generation is provided:
my $par2_prior=1-$par1_prior;
if($prior_admix==0){
    print "Running HMM "; print "ancestry_hmm -a 2 $par1_prior $par2_prior -p 0 -100 $par1_prior -p 1 -100 $par2_prior -s current.samples.list -i all.indivs.hmm.combined.filtered_focalchroms\n";
    system("ancestry_hmm -a 2 $par1_prior $par2_prior -p 0 -100 $par1_prior -p 1 -100 $par2_prior -s current.samples.list -i all.indivs.hmm.combined.filtered_focalchroms");
} else{

    if($par1_prior ge $par2_prior){
	print "Running HMM "; print "ancestry_hmm -a 2 $par1_prior $par2_prior -p 0 -100 $par1_prior -p 1 $initial_admix $par2_prior -s current.samples.list -i all.indivs.hmm.combined.filtered_focalchroms","\n";
    system("ancestry_hmm -a 2 $par1_prior $par2_prior -p 0 -100 $par1_prior -p 1 $initial_admix $par2_prior -s current.samples.list -i all.indivs.hmm.combined.filtered_focalchroms");
    } else{
	print "Running HMM "; print "ancestry_hmm -a 2 $par1_prior $par2_prior -p 0 $initial_admix $par1_prior -p 1 -100 $par2_prior -s current.samples.list -i all.indivs.hmm.combined.filtered_focalchroms","\n";
    system("ancestry_hmm -a 2 $par1_prior $par2_prior -p 0 $initial_admix $par1_prior -p 1 -100 $par2_prior -s current.samples.list -i all.indivs.hmm.combined.filtered_focalchroms");
    }#format properly depending on which parent is the major parent

}#if no gen is provided
