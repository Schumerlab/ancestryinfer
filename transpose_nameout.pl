#!perl
#
# 'transpose' swaps rows and columns in the given tab-delimited table.
#this version gets the names from an MSG ancestry.tsv file
#
my $infile=shift(@ARGV); chomp $infile;
open MYFILE, $infile or die "wrong format for infile\n";
   # $line=<MYFILE>;

my $outfile="$infile"."_transposed";
open OUT, ">$outfile";

while (<MYFILE>) {
    chomp;
    @line = split /\t/;
    $oldlastcol = $lastcol;
    $lastcol = $#line if $#line > $lastcol;
    for (my $i=$oldlastcol; $i < $lastcol; $i++) {
	$outline[$i] = "\t" x $oldlastcol;
    }
    for (my $i=0; $i <=$lastcol; $i++) {
    $outline[$i] .= "$line[$i]\t"
    }
}
for (my $i=0; $i <= $lastcol; $i++) {
    $outline[$i] =~ s/\s*$//g;
    print OUT $outline[$i]."\n";
}

#system("cut -f 1 loci > loci_temp");
#system("sed 1d loci_temp > loci_names");
#system("rm loci_temp loci_temp2");



