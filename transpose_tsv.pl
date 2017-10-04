#!perl

#MS modifying from: http://www1.cuni.cz/~obo/textutils/transpose

# 'transpose' swaps rows and columns in the given tab-delimited table.
#this version gets the names from an MSG ancestry.tsv file
#
my $infile=shift(@ARGV); chomp $infile;
open MYFILE, $infile or die "wrong format for infile\n";
   # $line=<MYFILE>;

my $outfile="$infile";
$outfile=~ s/_transposed//g;
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




