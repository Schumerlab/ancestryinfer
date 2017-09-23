#perl! -w

my $f1=shift(@ARGV); chomp $f1;
my $f2=shift(@ARGV); chomp $f2;
my $outfile=shift(@ARGV); chomp $outfile;

open HMM, ">$outfile";

open(F2,$f2); while (<F2>) { s/\r?\n//; @F=split /\t/, $_; $line2{$F[$col2]} .= "$_\n" }; $count2 = $.; open(F1,$f1); while (<F1>) { s/\r?\n//; @F=split /\t/, $_; $x = $line2{$F[$col1]}; if ($x) { $num_changes = ($x =~ s/^/$_\t/gm); print HMM $x; $merged += $num_changes } }
