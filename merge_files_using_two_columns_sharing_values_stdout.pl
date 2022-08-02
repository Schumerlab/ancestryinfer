#perl! -w

#chooses and merges rows in file2 based on a list in file1, using the values of two specified columns in each file
#useful for example for selecting particular sites from bed or vcf files where you need to match *both* the chromosome name and the site
#e.g.: ScyDAA6_2-HRSCAF_26     8585

#modified from the single file version on the FAS scriptome

if(@ARGV<6){
    print "perl merge_files_using_two_columns_sharing_values.pl file1_list column1_file1 column2_file1 file2_target column1_file2 column2_file2\n"; exit;
}#print usage

my $f1=shift(@ARGV); chomp $f1;
my $col1_1=shift(@ARGV); chomp $col1_1;
my $col1_2=shift(@ARGV); chomp $col1_2;

my $f2=shift(@ARGV); chomp $f2;
my $col2_1=shift(@ARGV); chomp $col2_1;
my $col2_2=shift(@ARGV); chomp $col2_2;

open(F2,$f2);
while (<F2>) {
    s/\r?\n//;
    @F=split /\t/, $_;
    $line2{$F[$col2_1]."_".$F[$col2_2]} .= "$_\n"
};
$count2 = $.;
open(F1,$f1);
while (<F1>) {
    s/\r?\n//;
    @F=split /\t/, $_;
    $x = $line2{$F[$col1_1]."_".$F[$col1_2]};
    if ($x) {
	$num_changes = ($x =~ s/^/$_\t/gm);
	print $x;
	$merged += $num_changes
    }
}


