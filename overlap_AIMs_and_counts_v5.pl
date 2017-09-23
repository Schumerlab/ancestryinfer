#perl! -w

my $aims=shift(@ARGV); chomp $aims;
open AIMS, $aims or die "cannot open AIMs file\n";

my $counts=shift(@ARGV); chomp $counts;
#!open COUNTS, $counts or die "cannot open counts file\n";

open OUT, ">"."$counts".".hmm";

my $match=""; my $combined="";
my $site=""; my $group="";
while(my $line =<AIMS>){

    chomp $line;

    my @data=split(/\t/,$line);
    $combined=$data[0]; chomp $combined;
    my @split_data=split(/_/,$combined);
    $group=$split_data[0]; chomp $group;
    $group =~ s/group//g;
    $site=$split_data[1]; chomp $site;

    my $ref_alt="$data[1]"."$data[2]"; chomp $ref_alt;
 
    $match=qx(grep -w $combined $counts); chomp $match;

    my $base1=""; my $base2=""; my $counts1=0; my $counts2=0;
    if(length($match)>0){

	my @matched_data=split(/\t/,$match);
	$base1=$matched_data[1]; chomp $base1;
	$base2=$matched_data[2]; chomp $base2;
	$counts1=$matched_data[3]; chomp $counts1;
	$counts2=$matched_data[4]; chomp $counts2;

	#!print "$base1\t$base2\t$counts1\t$counts2\n";

	my $order1="$base1"."$base2"; chomp $order1;
	my $order2="$base2"."$base1"; chomp $order2;

	if(($order1 == $ref_alt) or ($order2 == $ref_alt)){
	    print OUT "$group"."_"."$site\t20\t0\t0\t20\t0.0006\t$counts1\t$counts2\n";
	}#if match
	else{
	#!    print "$group\t$site\t1\t0\t0\t1\t0.0006\t0\t0\n";
	#!print "$combined\t$data[1]\t$data[2]\t0\t0\n";
	}#if some discordance, non-biallelic/AIM

    }#if there is a hits
    elsif(length($match)==0){
	    print OUT "$group"."_"."$site\t20\t0\t0\t20\t0.0006\t0\t0\n";
    }#if the AIM is not sampled in this individual

    $match="";

}#for all lines in AIMs file
