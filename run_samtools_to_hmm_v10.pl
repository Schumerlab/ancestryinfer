#perl! -w

if(@ARGV<9){

    print "perl run_samtools_to_hmm_v8.pl id_list genome1 genome2 read_length save_files max_align focal_chroms_file rate program_path quality\n"; exit;

}#usage

my $infile1=shift(@ARGV); chomp $infile1; my $trimmed="$infile1"."_trim";
system("cat $infile1 | perl -p -e 's/ +/\t/g' | cut -f 1 > $trimmed");
open IN, $trimmed or die "cannot open indiv id list\n";

my $genome1=shift(@ARGV); chomp $genome1;

my $genome2=shift(@ARGV); chomp $genome2;

my $read_length=shift(@ARGV); chomp $read_length;

my $save_files=shift(@ARGV); chomp $save_files;

my $max_align=shift(@ARGV); chomp $max_align;

my $focal_chroms=shift(@ARGV); chomp $focal_chroms;

my $rate=shift(@ARGV); chomp $rate; 

my $path=shift(@ARGV); chomp $path;

my $quality=shift(@ARGV); chomp $quality;

my $chrom_string="";
if($focal_chroms ne 0){
open FOCAL, $focal_chroms or die "cannot open focal_chroms_file\n";
while(my $tmp=<FOCAL>){
    chomp $tmp;
    if(length($chrom_string) eq 0){
    $chrom_string=$tmp;
    }#correctly generate names
    else{
    $chrom_string="$chrom_string".","."$tmp";
    }#append

}#generate string
} else{
    $chrom_string="allchrs";
}#only if defined
#print "$chrom_string\n";

open OUT, ">$infile1"."_bamlist";
my $outfile="$infile1"."_bamlist";

my $aims_list="current_aims_file"; chomp $aims_list;
open AIMSLIST, $aims_list or die "cannot open AIMs list\n";

my $aims=<AIMSLIST>; chomp $aims;
print "using AIMs in file $aims\n";

my @parfilelist=();
my @indivfilelist=();
while (my $id = <IN>){

    chomp $id;

    my $append="$chrom_string"; $append=~ s/,/_/g;
    print "file tag is $append";
    my $line1="$id"."_"."$append".".par1.sam";
    my $line2="$id"."_"."$append".".par2.sam";

    my $bam1="$line1".".bam";
    my $bam1_inter="$line1"."*".".bam";

    print "processing $bam1\nremoving if already exits\n";
    system("rm $bam1_inter");

    system("samtools fixmate -O bam $line1 $bam1");
    
    $sorted1="$bam1".".sorted";
    $sorted1=~ s/\.bam//g;

    $sorted1="$sorted1".".bam";
    
    print "sorting $sorted1\n";

    system("samtools sort $bam1 -o $sorted1");

    system("samtools index $sorted1");

    my $unique1 = "$sorted1".".unique.bam";
    $unique1 =~ s/sorted.bam.unique/sorted.unique/g;

    system("samtools view -b -q $quality $sorted1 > $unique1"); # this is for mapped reads with poor mapping quality

    print OUT "$unique1\n";

###now parent2

    my $bam2="$line2".".bam";
    my $bam2_inter="$line2"."*".".bam";

    print "processing $bam2\nremoving if already exits\n";
    system("rm $bam2_inter");

    system("samtools fixmate -O bam $line2 $bam2");

    $sorted2="$bam2".".sorted";
    $sorted2=~ s/\.bam//g;

    $sorted2="$sorted2".".bam";

    print "$sorted2\n";

    system("samtools sort $bam2 -o $sorted2");

    system("samtools index $sorted2");

    my $unique2 = "$sorted2".".unique.bam";
    $unique2 =~ s/sorted.bam.unique/sorted.unique/g;

    system("samtools view -b -q $quality $sorted2 > $unique2"); # this is for mapped reads with poor mapping quality                  

    print OUT "$unique2\n";

#####JOINT FILTERING
    my $par1_pass="$unique1"."_par1_passlist";
    my $par2_pass="$unique2"."_par2_passlist";

    $par1_pass=~ s/\//_/g;
    $par2_pass=~ s/\//_/g;

    $par1_pass=~ s/_read_1.fastq.gz.par1.sam.sorted.unique.bam//g;
    $par2_pass=~ s/_read_1.fastq.gz.par2.sam.sorted.unique.bam//g;

    system("samtools view -F 4 $unique1 | cut -f 1 > $par1_pass");
    system("samtools view -F 4 $unique2 | cut -f 1 > $par2_pass");

    my $pass_both="$par1_pass"."_both";
    $pass_both =~ s/_par1//g;

    ###intersect
    my $file1 = $par1_pass;
    my $file2 = $par2_pass;
    open F2, $file2 or die $!;
    open JOINT, ">$pass_both";
    while (<F2>) { $h2{$_}++ };
    open F1, $file1 or die;
    $total=$.; $printed=0;
    while (<F1>) { $total++; if ($h2{$_}) { print JOINT $_; $h2{$_} = ""; $printed++; } }

    ###filter
    my $finalbam1="$unique1";
    $finalbam1=~ s/sorted.unique/sorted.pass.unique/g;
    
    my $finalbam2="$unique2";
    $finalbam2=~s/sorted.unique/sorted.pass.unique/g;

    ###Adding max read filter here
    if($max_align > 0){
    print "limiting to $max_align alignments\n";
    my $tmp_list="$pass_both".".tmp";
    system("shuf -n $max_align $pass_both > $tmp_list");
    system("mv $tmp_list $pass_both");
    }#then subsample, otherwise leave pass list as is
    
    system("ngsutilsj bam-filter --whitelist $pass_both $unique1 $finalbam1");
    system("ngsutilsj bam-filter --whitelist $pass_both $unique2 $finalbam2");
    #system("bamutils filter $unique1 $finalbam1 -whitelist $pass_both");
    #system("bamutils filter $unique2 $finalbam2 -whitelist $pass_both");
   
    system("samtools index $finalbam1");
    system("samtools index $finalbam2");

    system("rm $par1_pass $par2_pass $pass_both");

#####VARIANT CALLING
    my $mpileup1 = "$unique1".".bcf";

    if($chrom_string ne 'allchrs'){
    system("bcftools mpileup -r $chrom_string -o $mpileup1 -f $genome1 $finalbam1");
    } else{
    system("bcftools mpileup -o $mpileup1 -f $genome1 $finalbam1");
    }#run on focal chromosomes or not

    my $vcf1 = "$unique1".".vcf.gz";

    system("bcftools call -mO z -o $vcf1 $mpileup1");

    system("gunzip $vcf1");

    $vcf1 = "$unique1".".vcf";

    system("perl $path/vcf_to_counts_non-colinear.pl $vcf1 $aims $path");

    my $counts1="$vcf1"."_counts";

    my $hmmsites1="$counts1".".hmm";

##modifying for bug here
##system("perl $path/overlap_AIMs_and_counts_v6.pl $aims $counts1 $rate");
    my $counts_present="$counts1"."_present";
    my $counts_absent="$counts1"."_absent";
    my $counts_combined="$counts1"."_combined";
    my $counts_combined_sorted="$counts1"."_combined_sorted";

    my $counts_bed="$counts1".".bed";
    system("cat $counts1 | perl -p -e 's/_/\t/g' | awk -v OFS=\'\\t\' \'\$1=\$1\"\\t\"\$2\' > $counts_bed");

    my $aims_bed="$aims".".mod.bed";

    system("bedtools intersect -a $aims_bed -b $counts_bed -wb -f 1 > $counts_present");
    system("bedtools intersect -a $aims_bed -b $counts_bed -wb -v -f 1 | awk -v OFS\=\'\\t\' \'\$6\=\$1\"\\t\"\$2\"\\t\"\$2\"\\t\"\$4\"\\t\"\$5\"\\t0\\t0\"\' > $counts_absent");

    system("cat $counts_present $counts_absent > $counts_combined");
    system("perl $path/merge_files_using_two_columns_sharing_values_stdout.pl $aims_bed 0 1 $counts_combined 0 1 | cut -f 1-5 --complement > $counts_combined_sorted");
    #MAKE SURE THE ORDER matches the aims file
    system("perl $path/reformat_counts_v2.pl $counts_combined_sorted $rate > $hmmsites1");

    my $hmm="$line1".".hmm.combined";
    $hmm =~ s/par1\.//g;

##split files before combining them
    my $append=$chrom_string;
    $append=~ s/,/_/g;
    my $parfilecurr="$hmm"."parental.format";
    my $indivfilecurr="$hmm".".pass.formatted";

    system("cut -f 1-6 $hmmsites1 > $parfilecurr");
    system("cut -f 7-8 $hmmsites1 > $indivfilecurr");

    push(@parfilelist, $parfilecurr);
    push(@indivfilelist,$indivfilecurr);

##cleanup intermediate files
    if($save_files==0){
	my $bai1="$sorted1".".bai"; my $bai2="$sorted2".".bai"; my $aimsvcf1="$vcf1".".aims"; $mod1="$vcf1".".mod"; 
    system("rm $line1 $line2 $bam1 $bam2 $sorted1 $sorted2 $unique1 $unique2 $finalbam1 $finalbam2 $mpileup1 $vcf1 $counts1 $bai1 $bai2 $aimsvcf1 $mod1 $counts_bed $counts_present $counts_absent $counts_combined");
    }#don't save files
}#done

my $append=$chrom_string;
$append=~ s/,/_/g;
my $listparfile="HMM.parental.files.list."."$infile1"."_"."$append"; $listparfile=~ s/\///g; $listparfile=~ s/\.\./\./g;
my $listhybfile="HMM.hybrid.files.list."."$infile1"."_"."$append"; $listhybfile=~ s/\///g; $listhybfile=~ s/\.\./\./g;
print "$listparfile\t$listhybfile\n";
open LISTPAR, ">"."$listparfile";
open LISTHYB, ">"."$listhybfile";

for my $j (0..scalar(@parfilelist)-1){

    print LISTPAR "$parfilelist[$j]\n";

}#parent file list

for my $l (0..scalar(@indivfilelist)-1){

    print LISTHYB "$indivfilelist[$l]\n";

}#parent file list 


