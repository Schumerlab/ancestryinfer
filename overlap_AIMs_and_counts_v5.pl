#perl! -w

###Snippet modifed from the FAS scriptome merge by columns snippet
my $f1=shift(@ARGV); chomp $f1;
my $f2=shift(@ARGV); chomp $f2;

open OUT, ">"."$f2".".hmm";

$col1=0;
$col2=0;

open(F2,$f2);
while (<F2>) {
    s/\r?\n//;
    @F=split /\t/, $_;
    $line2{$F[$col2]} .= "$_\n"
};
$count2 = $.;
open(F1,$f1);
while (<F1>) {
    s/\r?\n//;
    @F=split /\t/, $_;
    $x = $line2{$F[$col1]}; chomp $x;
    if ($x) {
	$x=~ s/\_[a-zA-Z]//g;
	my @elements_match=split(/\t/,$x);
       
	print OUT "$elements_match[0]"."\t"."20\t0\t0\t20\t0.00004\t"."$elements_match[1]"."\t"."$elements_match[2]\n";
	$merged += $num_changes
    }
    else{
	
	$_ =~ s/\_[a-zA-Z]//g;
	#print $_,"\t20\t0"."\n";
	print OUT "$_"."\t"."20\t0\t0\t20\t0.00004\t0\t0\n"
    }

}
