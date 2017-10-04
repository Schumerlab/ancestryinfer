#perl! -w

my $config=shift(@ARGV); chomp $config;
open CONFIG, $config or die "cannot open configuration file\n";

my $genome1=""; my $genome2=""; my $read_type=""; my $read_list=""; my $read_length=""; my $number_indiv_per_job="";
my $num_jobs=""; my $job1_submit=""; my $job2_submit=""; my $job3_submit=""; my $minor_prior=""; my $prefix="";
while (my $line=<CONFIG>){

    chomp $line;
    my @elements=split(/\=/,$line);
    
    if($line =~ /genome1/g){
	$genome1=$elements[1]; chomp $genome1;
	print "parent genome 1 is $genome1\n";
    }#define genome1
    if($line =~ /genome2/g){
	$genome2=$elements[1]; chomp $genome2;
	print "parent genome 2 is $genome2\n";
    }#define genome2
    if($line =~ /read_type/g){
	$read_type=$elements[1]; chomp $read_type;
	if(($read_type ne 'SE') && ($read_type ne 'PE')){
	    die "read type must be SE or PE\n";
	}
	print "read type is $read_type\n";
    }#define read type
    if($line =~ /read_list/g){
	$read_list=$elements[1]; chomp $read_list;
	print "read list is $read_list\n";
    }#define read  list
    if($line=~ /prop_genome_minor_parent/g){
	$minor_prior=$elements[1]; chomp $minor_prior;
	print "prior proportion of the genome from minor parent is $minor_prior\n";
    }#define priors
    if($line =~ /number_indiv_per_job/g){
	$number_indiv_per_job=$elements[1]; chomp $number_indiv_per_job;
       
	print "task being split into $number_indiv_per_job per job\n";
	$prefix="$read_list.";
	system("rm $prefix*");
	system("split -d -e -l $number_indiv_per_job $read_list $prefix"); 
	system("ls $prefix* > split_jobs_list");
    }#batch parameters
    if($line =~ /slurm_command_map/g){
	my @job1=split(/\#/,$line);
	$job1_submit="#"."$job1[1]"."\n"."#"."$job1[2]"."\n"."#"."$job1[3]"."\n"."#"."$job1[4]"."\n"."#"."$job1[5]"."\n"; chomp $job1_submit;
	print "cluster command for job1 is: $job1_submit\n";
    }#mapping command
    if($line =~ /slurm_command_variant_call/g){
        my @job2=split(/\#/,$line);
	$job2_submit="#"."$job2[1]"."\n"."#"."$job2[2]"."\n"."#"."$job2[3]"."\n"."#"."$job2[4]"."\n"."#"."$job2[5]"."\n"; chomp $job2_submit;
        print "cluster command for job2 is: $job2_submit\n";
    }#variant calling command 
    if($line =~ /slurm_command_hmm/g){
        my @job3=split(/\#/,$line);
	$job3_submit="#"."$job3[1]"."\n"."#"."$job3[2]"."\n"."#"."$job3[3]"."\n"."#"."$job3[4]"."\n"."#"."$job3[5]"."\n"; chomp $job3_submit;
        print "cluster command for job3 is: $job3_submit\n";
    }#mapping command 

}#read in the configuration file

####check if files exist
if(! -f $genome1){
    die "cannot find file $genome1\n";
}#genome2
if(! -f $genome2){
    die "cannot find file $genome2\n";
}#genome2
if(! -f $read_list){
    die "cannot find file $read_list\n";
}#read list

#####first make sure reference is indexed
my $g1bwt="$genome1".".sa";
my $g2bwt="$genome2".".sa";
if(! -f $g1bwt){
    system("bwa index $genome1");
}#index files for genome1 missing, build                                                                               
if(! -f $g2bwt){
    system("bwa index $genome2");
}#index files for genome2 missing, build

#####jobs array
my @jobs=();
open JOBS, "split_jobs_list";
while(my $tmp=<JOBS>){
    chomp $tmp; push(@jobs,$tmp);
}#collect array of jobs

####run analysis
##two types of cromwell files, one for SE and one for PE
my @slurm_ids_map=();
my $slurm_sam_string="";
my $hyb_string=""; my $par_string="";


    for my $j (0..scalar(@jobs)-1){
    my $current_job=$jobs[$j];
    #print "$current_job\n";
    open MAPSCRIPT, ">map_batch.sh";
    print MAPSCRIPT "$job1_submit\n";
    print MAPSCRIPT "perl run_map_v3.pl $current_job $genome1 $genome2 $read_type\n";
    
    my $id=qx(sbatch map_batch.sh);
    $id=~ s/Submitted batch job //g; chomp $id;
    push(@slurm_ids_map,$id);
    print "submitting mapping batch id $id\n";
    }#all mapping

    for my $m (0..scalar(@jobs)-1){
	my $current_job=$jobs[$m];

	my $hybfile="HMM.hybrid.files.list"."$current_job";
	my $parfile="HMM.parental.files.list."."$current_job";
	#save HMM file names for later
	if($m==0){
	    $hyb_string="$hybfile";
	    $par_string="$parfile";
	} else{
	    $hyb_string="$hyb_string"." "."$hybfile";
	    $par_string="$par_string"." "."$parfile"." ";
	}#generate file strings for next step

	open VARSCRIPT, ">samtools_batch.sh";
	print VARSCRIPT "$job2_submit\n";
	print VARSCRIPT "perl run_samtools_to_hmm_v5.pl $current_job $genome1 $genome2 150\n";

	my $map_depend=$slurm_ids_map[$m];

	my $id=qx(sbatch --dependency=afterok:$map_depend samtools_batch.sh);
	$id=~ s/Submitted batch job //g; chomp $id;
	if($m==0){
	    $slurm_sam_string="$id";
	} else{
	$slurm_sam_string="$slurm_sam_string".","."$id";
	}
	print "submitting variant batch id $id\n";

    }#all variants

    open HMMSCRIPT, ">hmm_batch.sh";
    print HMMSCRIPT "$job3_submit\n";
    print HMMSCRIPT "cat $hyb_string > HMM.hybrid.files.list\n";
    print HMMSCRIPT "cat $par_string >HMM.parental.files.list\n";
    print HMMSCRIPT "perl combine_all_individuals_hmm_v4.pl HMM.parental.files.list HMM.hybrid.files.list $minor_prior\n";
    print HMMSCRIPT "perl convert_rchmm_to_ancestry_tsv_v2.pl current.samples.list $read_list\n";
    print HMMSCRIPT "perl transpose_tsv.pl ancestry-probs-par1_transposed.tsv\n";
    print HMMSCRIPT "perl transpose_tsv.pl ancestry-probs-par2_transposed.tsv\n";
    system("sbatch --dependency=afterok:$slurm_sam_string hmm_batch.sh");
