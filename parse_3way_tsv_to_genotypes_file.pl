#perl! -w

if(@ARGV<2){
    print "perl parse_3way_tsv_ancestry.pl ancestry-file-stem posterior_prob_threshold > outfile\n"; exit;
}
#example stems:
#allchrs.tsv
#ScyDAA6-1508-HRSCAF-1794.tsv

my $stem=shift(@ARGV); chomp $stem;

my $thresh=shift(@ARGV); chomp $thresh;

my $par1="ancestry-probs-par1_"."$stem";
my $par2="ancestry-probs-par2_"."$stem";
my $par3="ancestry-probs-par3_"."$stem";
my $par1par2="ancestry-probs-par1par2_"."$stem";
my $par1par3="ancestry-probs-par1par3_"."$stem";
my $par2par3="ancestry-probs-par2par3_"."$stem";

open IN1, $par1 or die "cannot open $par1\n";
open IN2, $par2 or die "cannot open $par2\n";
open IN3, $par3 or die "cannot open $par3\n";
open IN4, $par1par2 or die "cannot open $par1par2\n";
open IN5, $par1par3 or die "cannot open $par1par3\n";
open IN6, $par2par3 or die "cannot open $par2par3\n";

my $header=<IN1>; my $junk2=<IN2>; my $junk3=<IN3>; my $junk4=<IN4>; my $junk5=<IN5>; my $junk6=<IN6>;
chomp $header;

$header="id\t$header";
$header=~ s/\t\t/\t/g;
#print "id\tpar1counts\tpar2counts\tpar3counts\tpar1\tpar2\tpar3\tpar1par2\tpar1par3\tpar2par3\n";

print "$header\n";

while((my $line1=<IN1>) && (my $line2=<IN2>) && (my $line3=<IN3>) && (my $line4=<IN4>) && (my $line5=<IN5>) && (my $line6=<IN6>)){

    chomp $line1; chomp $line2; chomp $line3; chomp $line4; chomp $line5; chomp $line6;

    my @elements1=split(/\t/,$line1);
    my @elements2=split(/\t/,$line2);
    my @elements3=split(/\t/,$line3);
    my @elements4=split(/\t/,$line4);
    my @elements5=split(/\t/,$line5);
    my @elements6=split(/\t/,$line6);

    my $id=$elements1[0]; chomp $id;
    my $par1=0; my $par2=0; my $par3=0; my $par1par2=0; my $par1par3=0; my $par2par3=0; my $total=0;
    my $geno="NA";
    print "$id\t";
    for my $i (1..scalar(@elements1)-1){

	my $focal1=$elements1[$i];
	my $focal2=$elements2[$i];
	my $focal3=$elements3[$i];
	my $focal4=$elements4[$i];
	my $focal5=$elements5[$i];
	my $focal6=$elements6[$i];

	#print "$focal1\t$focal2\t$focal3\n";

	if($focal1 >= $thresh){$par1++; $total++; $geno=0;}
	if($focal2 >= $thresh){$par2++; $total++; $geno=2;}
	if($focal3 >= $thresh){$par3++; $total++; $geno=5;}
	if($focal4 >= $thresh){$par1=$par1+0.5; $par2=$par2+0.5; $par1par2++; $total++; $geno=1;}
	if($focal5 >= $thresh){$par1=$par1+0.5; $par3=$par3+0.5; $par1par3++; $total++; $geno=3;}
	if($focal6 >= $thresh){$par2=$par2+0.5; $par3=$par3+0.5; $par2par3++; $total++; $geno=4;}

	if($i<(scalar(@elements1)-1)){
	    print "$geno\t";
	} else{
	    print "$geno\n";
	}

    }#for all elements

}#for all individuals
