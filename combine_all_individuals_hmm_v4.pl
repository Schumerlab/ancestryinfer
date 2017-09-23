#perl! -w

if(@ARGV<2){
    print "perl combine_all_individuals_hmm.pl parental_list individual_list\n";
}

my $parental=shift(@ARGV); chomp $parental;
open PAR, $parental or die "cannot open parental list file\n";

my $list=shift(@ARGV); chomp $list;
open LIST, $list or die "cannot open individual list file\n";

my $par1_prior=shift(@ARGV); chomp $par1_prior;
if(length($par1_prior)==0){ die "parental priors required\n";}

open OUT, ">current.samples.list";

my $counter=0;
my $string="";
while ((my $line1=<LIST>) && (my $line2=<PAR>)){

    $counter++;
    chomp $line1; chomp $line2;

    if($counter==1){
       $string="$line2"." "."$line1"." ";
    } else{
	$string="$string"." "."$line1"." ";
    }

    my $posterior="$line1".".posterior";
    print OUT "$posterior\t2\n";

}

print "$string\n";
system("paste $string > all.indivs.hmm.combined");

open CHMM, "all.indivs.hmm.combined" or die "cannot open combined HMM input data\n";
open FILTER, ">all.indivs.hmm.combined.filtered";
my $total=0;
while(my $line=<CHMM>){

    chomp $line;
    my @dataelements=split(/\t/,$line);
    $total=0;
    for my $i (7..scalar(@dataelements)-1){
        my $focal=$dataelements[$i]; chomp $focal;
        $total=$total+$focal;
    }

    if($total >0){
        print FILTER "$line\n";
    }

}

system("perl -pi -e 's/\\//_/g' current.samples.list"); 

my $par2_prior=1-$par1_prior;
system("ancestry_hmm -a 2 $par2_prior $par1_prior -p 0 -100 -$par2_prior -p 1 -100 -$par1_prior -s current.samples.list -i all.indivs.hmm.combined.filtered") 
