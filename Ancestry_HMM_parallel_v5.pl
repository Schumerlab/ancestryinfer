#perl! -w

my $config=shift(@ARGV); chomp $config;
open CONFIG, $config or die "cannot open configuration file\n";

my $genome1=""; my $genome2=""; my $read_type=""; my $read_list=""; my $read_length=""; my $number_indiv_per_job=""; my $initial_admix=0; my $focal_chrom=0; my $pp=0.9;
my $num_jobs=""; my $job1_submit=""; my $job2_submit=""; my $job3_submit=""; my $minor_prior=""; my $prefix=""; my $error_prior=""; my $max_align=""; my $rec_M_per_bp=0.00000003; my $path=""; my $job_submit="";
my $provide_AIMs=""; my $provide_counts=""; my $parental_counts_status=0; my $save_files=0; my $quality=30;
while (my $line=<CONFIG>){

    chomp $line;
    my @elements=split(/\=/,$line);
    
    if($line =~ /genome1=/g){
	$genome1=$elements[1]; chomp $genome1;
	print "parent genome 1 is $genome1\n";

	if($genome1 =~ /\.gz/g){
	    print "unzipping $genome1\n";
	    system("gunzip $genome1");
	    $genome1 =~ s/\.gz//g;
	}#unzip and rename

    }#define genome1
    if($line =~ /genome2=/g){
	$genome2=$elements[1]; chomp $genome2;
	print "parent genome 2 is $genome2\n";

	if($genome2 =~ /\.gz/g){
            print "unzipping $genome2\n";
            system("gunzip $genome2");
	    $genome2  =~ s/\.gz//g;
        }#unzip and rename

    }#define genome2
    if($line =~ /program_path=/g){
	$path=$elements[1]; chomp $path;
    }#path
    if($line =~ /read_type/g){
	$read_type=$elements[1]; chomp $read_type;
	if(($read_type ne 'SE') && ($read_type ne 'PE')){
	    die "read type must be SE or PE\n";
	}
	print "read type is $read_type\n";
    }#define read type
    if($line =~ /read_length/g){
	$read_length=$elements[1]; chomp $read_length;
	print "expected read length is $read_length\n";
    }#save read length
    if($line =~ /mapping_quality/g){
	my $tmpquality=$elements[1]; chomp $tmpquality;
	if(length($tmpquality)>0){
	    $quality=$tmpquality;
	}#reset quality
	print "using mapping quality cutoff of $quality\n";
    }#quality cutoff
    if($line =~ /read_list/g){
	$read_list=$elements[1]; chomp $read_list;
	print "read list is $read_list\n";
    }#define read  list
    if($line=~ /prop_genome_genome1_parent/g){
	$minor_prior=$elements[1]; chomp $minor_prior;
	print "prior proportion of the genome from parent 1 is $minor_prior\n";
    }#define priors
    if($line=~ /gen_initial_admix/g){
	$initial_admix=$elements[1]; chomp $initial_admix;
	if(length($initial_admix)==0){$initial_admix=0;} else{print "prior time of initial admixture is $initial_admix\n";}
    }#define prior admixture expectation
    if($line=~ /per_site_error/g){
	$error_prior=$elements[1]; chomp $error_prior;
	if(length($error_prior)==0){$error_prior=0;}
    }#define the expected error parameter for the HMM
    if($line=~ /max_alignments/g){
	$max_align=$elements[1]; chomp $max_align;
	if(length($max_align)==0){$max_align=0;}
    }#set the number of alignments to retain in the sam file
    if($line=~ /focal_chrom_list/g){
	$focal_chrom=$elements[1]; chomp $focal_chrom;
	if(length($focal_chrom)==0){$focal_chrom=0;}
    }#define focal chromosomes for running HMM
    if($line=~ /retain_intermediate_files/g){
	$save_files=$elements[1]; chomp $save_files;
    }#define save files
    if($line=~ /job_submit_command/g){
	$job_submit=$elements[1]; chomp $job_submit;
    }#define submit command as sbatch or bash
    if($line=~ /posterior_thresh/g){
	$pp=$elements[1]; chomp $pp;
	if(length($pp) eq 0){$pp=0.9}
    }#read in pp thresh for rec events
    if($line=~ /rec_M_per_bp/g){
	$rec_M_per_bp=$elements[1]; chomp $rec_M_per_bp;
    }#set per bp recombination rate
    if($line =~ /number_indiv_per_job/g){
	$number_indiv_per_job=$elements[1]; chomp $number_indiv_per_job;
       
	print "task being split into $number_indiv_per_job per job\n";

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
    if($line =~ /provide_AIMs/){
	$provide_AIMs=$elements[1]; chomp $provide_AIMs;
	
	if($provide_AIMs =~ /\.gz/g){
            print "unzipping $provide_AIMs\n";
	    system("gunzip $provide_AIMs");
	    $provide_AIMs =~ s/\.gz//g;
	}#unzip and rename

    }#aims list if provided
    if($line =~ /provide_counts/){
        $provide_counts=$elements[1]; chomp $provide_counts;
	if(length($provide_counts)>0){
	$parental_counts_status=$provide_counts;
	} else{
	$parental_counts_status=0;
	}#id variable

	if($provide_counts =~ /\.gz/g){
	    print "unzipping $provide_counts\n";
            system("gunzip $provide_counts");
            $provide_counts =~ s/\.gz//g;
	}#unzip and rename

    }#aims list if provided 

}#read in the configuration file

####generate job split
if($job_submit eq 'sbatch'){
    $prefix="$read_list.";
    system("rm $prefix*"); #remove anything from past runs
    system("split -d -e -l $number_indiv_per_job $read_list $prefix");
    system("ls $prefix* > split_jobs_list");
} else{
    system("ls $read_list > split_jobs_list");
}#split or don't

####check if files exist
if(length($path) < 1){die "Please set a path to the Ancestry_HMM_pipeline folder in the configuration file\n";} #path check 
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


#####generate AIMs list or identify use-defined AIMs list
my $aims="all_AIMs_"."$genome1"."_"."$genome2";
$aims=~ s/\.\///g;
$aims=~ s/\//_/g;

if(length($provide_AIMs)==0){
    system("perl $path/identify_AIMs_two_genomes.pl $genome1 $genome2 > $aims");
    open AIMSFILE, ">current_aims_file";
    print AIMSFILE "$aims\n";
}#if aims file and key do not exist, write them
elsif(-f $provide_AIMs){

    system("cp $provide_AIMs $aims");
    system("perl -pi -e 's/\t/_/g' $aims"); #reformat for downstream compatibility
    open AIMSFILE, ">current_aims_file";
    print AIMSFILE "$aims\n";

}#aims file defined and exists
elsif(! -f $provide_AIMs){

    die "cannot open the user defined AIMs file $provide_AIMs\n";

}#aims file defined but does not exist

#####reformat AIMs file for easy overlap
my $aims_bed="$aims".".mod";
system("cat $aims | perl -p -e 's/_/\t/g' | awk -v OFS=\'\\t\' \'\$1=\$1\"\_\"\$2\' > $aims_bed");

#####jobs array
my @jobs=();
open JOBS, "split_jobs_list";
while(my $tmp=<JOBS>){
    chomp $tmp; push(@jobs,$tmp);
}#collect array of jobs

##generate run-specific file names                                                                                      
my $tag="";
if($focal_chrom ne 0){
    open CHR, "$focal_chrom";
    while (my $chrnames=<CHR>){
        chomp $chrnames;
	$tag="$tag"."_"."$chrnames";
    }#generate chrom specific tags                                                                                      
} else{
    $tag="_allchrs";
}#name depends on focal_chrom status

####run analysis
##two types of approaches, one for SE and one for PE
my @slurm_ids_map=();
my $slurm_sam_string="";
my $hyb_string=""; my $par_string="";


    for my $j (0..scalar(@jobs)-1){
    my $current_job=$jobs[$j];
    #print "$current_job\n";
    my $mapscript="map_batch"."$j".".sh";
    open MAPSCRIPT, ">$mapscript";
    print MAPSCRIPT "$job1_submit\n";
    print MAPSCRIPT "perl $path/run_map_v3.pl $current_job $genome1 $genome2 $read_type $tag\n";

    if($job_submit eq 'sbatch'){
    my $id_temp=qx(sbatch $mapscript); chomp $id_temp;
    my @idarray=split(/\n/,$id_temp);
    $id=$idarray[0]; chomp $id;
    $id=~ s/\D//g;
    push(@slurm_ids_map,$id);
    print "submitting mapping batch id $id\n";
    } else{
	system("bash $mapscript");
    }#submit parallel or sequential

    }#all mapping



    for my $m (0..scalar(@jobs)-1){
	my $current_job=$jobs[$m];

	my $hybfile="HMM.hybrid.files.list."."$current_job"."$tag"; $hybfile=~ s/\///g; $hybfile =~ s/\.\./\./g;  
	my $parfile="HMM.parental.files.list."."$current_job"."$tag"; $parfile=~ s/\///g; $parfile =~ s/\.\./\./g; 

	print "$hybfile\t$parfile\n";

	#save HMM file names for later
	if($m==0){
	    $hyb_string="$hybfile";
	    $par_string="$parfile";
	} else{
	    $hyb_string="$hyb_string"." "."$hybfile";
	    $par_string="$par_string"." "."$parfile"." ";
	}#generate file strings for next step

	my $samscript="samtools_batch"."$m".".sh";
	open VARSCRIPT, ">$samscript";
	print VARSCRIPT "$job2_submit\n";
	print VARSCRIPT "perl $path/run_samtools_to_hmm_v8.pl $current_job $genome1 $genome2 $read_length $save_files $max_align $focal_chrom $rec_M_per_bp $path $quality\n";

	my $map_depend=$slurm_ids_map[$m]; 

	if($job_submit eq 'sbatch'){
	my $id_temp=qx(sbatch --dependency=afterok:$map_depend $samscript); chomp $id_temp;
	my @idarray=split(/\n/,$id_temp);
	$id=$idarray[0]; chomp $id;
	$id=~ s/\D//g;
	if($m==0){
	    $slurm_sam_string="$id";
	} else{
	$slurm_sam_string="$slurm_sam_string".","."$id";
	}
	print "submitting variant batch id $id\n";
	} else{
	    system("bash $samscript");
	}#submit parallel or sequential 

    }#all variants

##output files
my $final_file1="ancestry-probs-par1_transposed"."$tag".".tsv"; my $final_file2="ancestry-probs-par2_transposed"."$tag".".tsv";

    print "output files appended with $tag\n";
    open HMMSCRIPT, ">hmm_batch.sh";
    print HMMSCRIPT "$job3_submit\n";
    print HMMSCRIPT "cat $hyb_string > HMM.hybrid.files.list"."$tag"."\n";
    print HMMSCRIPT "cat $par_string >HMM.parental.files.list"."$tag"."\n";
    print HMMSCRIPT "perl $path/combine_all_individuals_hmm_v5.pl HMM.parental.files.list"."$tag HMM.hybrid.files.list"."$tag $minor_prior $parental_counts_status $initial_admix $focal_chrom $read_length $error_prior $tag\n";
#print "perl combine_all_individuals_hmm_v5.pl HMM.parental.files.list HMM.hybrid.files.list $minor_prior $parental_counts_status $initial_admix $focal_chrom $read_length $error_prior $tag\n";
    print HMMSCRIPT "perl $path/convert_rchmm_to_ancestry_tsv_v3.pl current.samples.list current.samples.read.list $save_files $focal_chrom\n";
    print HMMSCRIPT "perl $path/transpose_tsv.pl $final_file1\n";
    print HMMSCRIPT "perl $path/transpose_tsv.pl $final_file2\n";
    print HMMSCRIPT "rm split_jobs_list\n"; #cleanup split read lists

$final_file1=~ s/_transposed//g; $final_file2=~ s/_transposed//g;
$rec_geno="$final_file1"."_rec.txt";
$rec_geno=~ s/-par1//g;
    print HMMSCRIPT "perl $path/parsetsv_to_genotypes_v2.pl $final_file1 $final_file2 $pp $rec_geno\n";
    print HMMSCRIPT "Rscript $path/identify_intervals_ancestryinfer.R $rec_geno $path\n";

    print HMMSCRIPT "rm ancestry-probs_"."*",".tsv_rec.txt\nrm ancestry-probs_"."*".".tsv_rec.txt_header\nrm ancestry-probs_"."*".".tsv_rec.txt_transposed"."*"."\n"; #cleanup intermediate files         


if($job_submit eq 'sbatch'){
    system("sbatch --dependency=afterok:$slurm_sam_string hmm_batch.sh");
} else{
    system("bash hmm_batch.sh");
}#submit parallel or sequential  
