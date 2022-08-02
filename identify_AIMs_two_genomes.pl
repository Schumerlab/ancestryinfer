#perl! -w

my $infile1=shift(@ARGV); chomp $infile1;
open IN1, $infile1 or die "cannot open fasta species1\n";

my $infile2=shift(@ARGV); chomp $infile2;
open IN2, $infile2 or die "cannot open fasta species2\n";

my $sites1="";
my $sites2="";
my $pos=0;
my $lg="";
my $focal1="";
my $focal2="";
while((my $line1=<IN1>) && (my $line2=<IN2>)){
    chomp $line1; chomp $line2;
    if($line1=~ />/){
	$lg=$line1;
	$lg =~ s/>//g;
	#print "$lg\n";
	$pos=0;
    } else{

    my @elements1=split(//,$line1);
    my @elements2=split(//,$line2);

    if(scalar(@elements1) != scalar(@elements2)){
	die;
    }

    
    for my $i (0..scalar(@elements1)-1){
	my @bparray=();

	$pos=$pos+1;

	$focal1=$elements1[$i]; chomp $focal1;
	$focal2=$elements2[$i]; chomp $focal2;

	push(@bparray, $focal1);
	push(@bparray, $focal2);
	
#!	my @sortedbp=sort { lc($a) cmp lc($b) } @bparray;

	if(($focal1 ne $focal2) && ($focal1 !~ /[RYSWKMN]/) && ($focal2 !~ /[RYSWKMN]/)){
	    print "$lg"."\t"."$pos"."\t"."$bparray[0]"."\t"."$bparray[1]\n";

	}#sites are AIMs in genome

    }#for all sites

    }#non-header line

}#for all lines
