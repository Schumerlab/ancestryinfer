#perl! -w

my $infile=shift(@ARGV); chomp $infile;
open IN, $infile or die "cannot open counts file\n";

my $rec=shift(@ARGV); chomp $rec;

my $prev=0; my $dist=0;
while(my $line=<IN>){

    chomp $line;
    my @elements=split(/\t/,$line);

    my $current=$elements[1]; $dist=$current-$prev; $prev=$current;
    #print "$dist\n";

    my $ref1=$elements[3]; chomp $ref1;
    my $ref2=$elements[8]; chomp $ref2;
    my $alt1=$elements[4]; chomp $alt1;
    my $alt2=$elements[9]; chomp $alt2;

    #print "$ref1\t$ref2\t$alt1\t$alt2\n";

    if(($ref1 eq $ref2) && ($alt1 eq $alt2)){
    print "$elements[0]"."_"."$elements[1]"."\t20\t0\t0\t20\t";
    my $rate=$dist*$rec;
    print "$rate\t"."$elements[10]"."\t"."$elements[11]"."\n";
    }
    elsif(($ref1 eq $alt2) && ($alt1 eq $ref2)){
    print "$elements[0]"."_"."$elements[1]"."\t20\t0\t0\t20\t";
    my $rate=$dist*$rec;
    print "$rate\t"."0"."\t"."0"."\n";
    }
    

}
