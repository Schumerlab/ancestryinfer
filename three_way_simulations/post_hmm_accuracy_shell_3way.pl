#perl! -w

if(@ARGV<4){
    print "perl post_hmm_accuracy_shell.pl tsv_stem_name bed_list posterior_thresh path_to_Simulator_scripts\n"; exit;
}

#example command:
#perl post_hmm_accuracy_shell_3way.pl allchrs.tsv bed_list 0.9 ./

#bed list format (auto output of the simulator script, e.g. admixture_SELAM_sim_bed_list ):
#bed_hap1\tbed_hap2

#tsv stem name format example:
#allchrs.tsv


my $stem=shift(@ARGV); chomp $stem;

my $bed_list=shift(@ARGV); chomp $bed_list;
open IN, $bed_list or die "cannot open bed list\n";

my $pp=shift(@ARGV); chomp $pp;

my $path=shift(@ARGV); chomp $path;

my $genos="genotypes"."_"."simulation_accuracy";
system("perl $path/parse_3way_tsv_to_genotypes_file.pl $stem $pp > $genos");

my $outfile="results_summary_"."$bed_list"; 
open OUT, ">results_summary_"."$bed_list";
while(my $line=<IN>){
    chomp $line;
    my @elements=split(/\t/,$line);
    my $bed1=$elements[0];
    my $bed2=$elements[1];

    my @id_data=split(/_tracts_/,$bed1);
    my $id=$id_data[0]; chomp $id;
    print "$id\n";

    my $recode="$bed1"."_combined_recode";
    system("Rscript recode_bed_files.R $bed1 $bed2");

    my $current=qx(Rscript $path/Determine_accuracy_3way.R $recode $id $genos $path); chomp $current;
    print OUT "$current\n";

}#for all individuals

system("rm accuracy_*genotypes_file*");



