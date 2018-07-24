#!/usr/bin/perl

use strict;
use lib ''; #put here the full path of the folder with the perl scripts
use Create_cmd;


#########################################################################
# Pipeline running of raw read files (fastq.gz) to frequency files
# Begins with an input direcory containing fastq.gz files
# 1. Convert fastq.gz to fasta (mask low Q bases - V3), and
# 2. Split all fasta files into N equally sized smaller fasta files, with ~50K reads in each
# 3. Run formatdb on each of N files above, and blast against ref seq
# 4. run base calling script on each blast file above (output-freq files)
# 5. Join output of all freq files above
#########################################################################


die "usage pipeline_runner.pl  <input directory with fastq.gz files> <output dir> <reference genome seq (fasta)> <start at step number> <end at step number> <type of input files, optional f if fastq and not zipped files> <refer to gaps? Y/N default Y> <NGS/Cirseq?  type 1 for NGS (min num repeats=1) and >1 for AccuNGS> <Q-score cutoff, default =30 for NGS> <% id for blast, default=85><E value for blast, default=1e-7>\n
1. Convert fastq.gz to fasta & split all fasta files into N equally sized smaller fasta files (50K reads per file)\n
2. Run formatdb on each of N files above, and blast against ref seq\n
3. run base calling script on each blast file above (output-freq files)\n
4. Join output of all freq files above

" unless (scalar(@ARGV)>=5);


my $current_version_path = '';  #put here the full path of the folder with the perl scripts
my $in_dir = $ARGV[0];

$in_dir.="/" unless ($in_dir =~ /\/$/);

my $out_dir = $ARGV[1];
$out_dir.="/" unless ($out_dir =~ /\/$/);

unless (-e $out_dir) {system("mkdir $out_dir");}


my $file_to_write_in_dir = $out_dir."fastq_dir";
open FA, ">$file_to_write_in_dir" or die "cannot open file $file_to_write_in_dir\n";
print FA "input dir for fastq files is $in_dir\n";
close FA;

my $ref_genome = $ARGV[2];
unless (-e $ref_genome) {die "error, reference genome $ref_genome does not exist\n"};

my $err_file="${out_dir}/errors.txt";

my $start_stage=$ARGV[3];

my $end_stage=$ARGV[4];

my $type_file="z"; # default - gzipped fastq files

if (defined $ARGV[5]) {
    $type_file=$ARGV[5];
}

my $do_gaps="Y"; # default - Y, refer to gaps.
if (defined $ARGV[6]) {
    $do_gaps=$ARGV[6];
}

my $min_num_repeats=$ARGV[7]; 
if ($min_num_repeats==1) {
    print "Running NGS mapping\n";
}
elsif ($min_num_repeats>1) {
    print "Running CirSeq mapping\n";

}
else {
    die "min number of repeats should be either 1 or bigger, it is now $min_num_repeats\n";
}



# default quality cutoff per base is 30
my $q_cutoff=30;

if (defined $ARGV[8]) {
    $q_cutoff=$ARGV[8];
}

my $pcID_blast = 85;


if (defined $ARGV[9]) {

    $pcID_blast=$ARGV[9];

}

my $evalue = 1e-7;


if (defined $ARGV[10]) {

    $evalue=$ARGV[10];

}



die "unexpected error, start stage $start_stage is larger than end stage $end_stage\n" if ($start_stage>$end_stage);


my $scripts_dir = $current_version_path;
my $ascii_file = "$scripts_dir/ascii_table_processed.txt";
my $blast_dir = ""; #put here the full path of the BLAST executables (or put them in the PATH variable)



my $number_of_reads_per_fasta_file=50000;
my $list_files=$out_dir."list.fastq.txt";
my $num_fasta_files_in_this_run=-1;


my $sleep_quantum=60;
my $sleep_max=1200000; 

open ERR, ">$err_file" or die "cannot open file $err_file\n";
&main;
close ERR;



sub main {
    for my $i($start_stage..$end_stage) {
	   my $num_file=&process($i);
    }
}

sub process {
    my $num=shift;

    &splitFastq if ($num==1);
    &blast_fasta if ($num==2);
    &base_call if ($num==3);    
    if ($num==4) {
		&join_files;
		&wrap_up;
    }
}



sub splitFastq {

    print "splitFastq\n";

    my $cmd="";
    if ($type_file eq "z") {
		$cmd="ls -l $in_dir"."*.gz \| awk \'{print \$NF}\' \| awk -F \"/\" \'{print \$NF}\' \| perl -lane \'s/\.fastq\.gz//\;print;' >$list_files";
    }
    else {
	   $cmd="ls -l $in_dir"."*.fastq \| awk \'{print \$NF}\' \| awk -F \"/\" \'{print \$NF}\' \| perl -lane \'s/\.fastq//\;print;' >$list_files";
    }
	system($cmd);

    my $num_files=&get_number_of_fasta_files_there_should_be($list_files);
    my $alias="toFastaAndSplit";
    my $cmd_file=$out_dir."splitFastq.cmd";

    my $cmd1="";



    if ($num_files==1) {	# if only one file is present, ARRAY_INDEX is not defined
        $cmd1 = "INFILE\=\$\(awk \"NR\=\=1\" $list_files\)\n";
    }
    else {
       $cmd1 = "INFILE\=\$\(awk \"NR\=\=\$PBS_ARRAY_INDEX\" $list_files\)\n";
    }

    print "fasta cmd1: $cmd1\n";

    my $cmd2="";
    if ($type_file eq "z") {
		$cmd2="perl $scripts_dir/toFastaAndSplit.pl $in_dir\$INFILE.fastq.gz $out_dir\$INFILE $number_of_reads_per_fasta_file\n";
    }
    else {
		$cmd2="perl $scripts_dir/toFastaAndSplit.pl $in_dir\$INFILE.fastq $out_dir\$INFILE $number_of_reads_per_fasta_file f\n";
	}

    print "cmd2 is $cmd2\n";

    Create_cmd::create_cmd_file($cmd_file,$alias,$num_files,4,join(" ",$cmd1,$cmd2));
    $cmd="qsub $cmd_file";

   	 

#Your job-array 8416810.1-4:1 ("toFastaAndSplit") has been submitted

    my $stdin=`$cmd`;
    &sleep_manager(1,$stdin,$alias);

    sleep(60); # to ensure the glob test works... need to wait
    my @fasta_files=glob($out_dir."*fasta");

    if ($num_files != scalar(@fasta_files)) {
	print ERR "error in toFastaAndSplit, number of available fasta files in directory $out_dir = ".scalar(@fasta_files).", should be $num_files\n ";
    }
    if (scalar(@fasta_files)==0) {
	   die "error in fastq2fasta, no fasta files created\n"; 
    } 
} ##splitFastq



sub blast_fasta {

    print "blast_fasta\n";
    my $list_parts_fasta_files = $out_dir."list_parts_fasta_files.txt";
	
    # all files with size >0
    my $glob_test=$out_dir."*part*fasta";
    my @test=glob($glob_test);

    die "Error, no files produced from splitfasta, glob tests is $glob_test\n" if (scalar(@test)==0);
    my $cmd="ls -lrt $glob_test \| awk \'\$5>0\' \| awk \'{print \$NF}\' \| awk -F \"/\" \'{print \$NF}\' >$list_parts_fasta_files";
    system($cmd);

    $cmd="wc $list_parts_fasta_files \| awk \'{print \$1}\'";

    my $num_files=`$cmd`;
    my $alias="blast";
    my $cmd_file=$out_dir."blastfasta.cmd";

    my $cmd1="";

    print "blast\n";

    if ($num_files==1) {
        $cmd1 = "INFILE\=\$\(awk \"NR\=\=1\" $list_parts_fasta_files\)\n";
    }
    else {
       $cmd1 = "INFILE\=\$\(awk \"NR\=\=\$PBS_ARRAY_INDEX\" $list_parts_fasta_files\)\n";
    }

    print "cmd1 blast: $cmd1\n";

    print "blast: list_parts_fasta_files: $list_parts_fasta_files\n";

    my $cmd2="$blast_dir/makeblastdb -in $out_dir\$INFILE -dbtype nucl\n";
	my $cmd3="$blast_dir/blastn -query $ref_genome -task blastn -db $out_dir\$INFILE -outfmt \"6 sseqid qstart qend qstrand sstart send sstrand length btop\" -num_alignments 1000000 -dust no -soft_masking F -perc_identity $pcID_blast -evalue $evalue -out $out_dir\$INFILE.blast";

    print "cmd3 blast: $cmd3\n";

    Create_cmd::create_cmd_file($cmd_file,$alias,$num_files,4,join(" ",$cmd1,$cmd2,$cmd3));
    $cmd="qsub $cmd_file";

    my $stdin=`$cmd`;
    
    &sleep_manager(3,$stdin,$alias);
    
    my @files=glob("$out_dir*.fasta");
    my @files_blast=glob("$out_dir*.blast");

    if (scalar(@files)!=scalar(@files_blast) ) {
	   print ERR "number of blast output files".scalar(@files_blast) ." not compatible with number of input fasta files".scalar(@files)."\n";
    }

    if (scalar @files_blast == 0) {
	   die "error in blastfasta: no blast output files produced\n";
    }
}# blast_fasta



sub base_call {

    print "base_call\n";

    my $list_blast_results_files = $out_dir."list_blast_results_files.txt";
	my $cmd="ls -lc $out_dir"."*blast \| awk \'\$5>0\' \| awk \'{print \$NF}\' | sort >$list_blast_results_files";
    system($cmd);

    my $list_qual_files = $out_dir."list_qual_files.txt";
    my $qcmd="ls -lc $out_dir"."*qual \| awk \'\$5>0\' \| awk \'{print \$NF}\' | sort >$list_qual_files";
    system($qcmd);

    $cmd="wc $list_blast_results_files \| awk \'{print \$1}\'";
    my $num_files=`$cmd`;
    my $alias="basecall";
    my $cmd_file=$out_dir."basecall.cmd";

    print "base_call: list_parts_fasta_files:".$list_blast_results_files;

    my $cmd1="";
    my $cmd2="";

    if ($num_files==1) {
        $cmd1="INFILE\=\$\(awk \"NR\=\=1\" $list_blast_results_files\)\n";
        $cmd2="FASTQ\=\$\(awk \"NR\=\=1\" $list_qual_files\)\n";
    }
    else {
       	$cmd1="INFILE\=\$\(awk \"NR\=\=\$PBS_ARRAY_INDEX\" $list_blast_results_files\)\n";
     	$cmd2="FASTQ\=\$\(awk \"NR\=\=\$PBS_ARRAY_INDEX\" $list_qual_files\)\n";
    }

    my $cmd3="perl $scripts_dir/base_call_and_freqs.pl \$INFILE \$FASTQ $ref_genome \$INFILE\.freqs $do_gaps $min_num_repeats $q_cutoff\n";
    print "cmd3 basecall: $cmd3\n";
    my $mem_request = 8; 
    Create_cmd::create_cmd_file($cmd_file,$alias,$num_files,$mem_request,join(" ",$cmd1,$cmd2,$cmd3));

    $cmd="qsub $cmd_file";

    my $stdin=`$cmd`;
    &sleep_manager(4,$stdin,$alias);

    my @files_blast=glob("$out_dir*.blast");
    my @files_freqs=glob("$out_dir*.freqs");

    if (scalar(@files_freqs)!=scalar(@files_blast) ) {
	   print ERR "number of blast output files ".scalar(@files_blast) ." not compatible with number of freqs files created: ".scalar(@files_freqs)."\n";
    }

    if (scalar @files_freqs == 0) {
	   die "error in blastfasta: no freqs output files produced\n";
    }
} # base_call





sub join_files {

    print "join_files\n";

    my $list_freqs_files = $out_dir."*part*.freqs";

	# list all output blast files with file size > 0
    my $cmd="ls -l $list_freqs_files \| awk \'{print \$NF}\' \| awk -F \"/\" \'{print \$NF}\' \| head -1 \| awk -F \"\.\" \'{print \$1}\'";

    my $prefix=`$cmd`;
    chomp $prefix;

    $prefix =~ s/\.fastq//;
    $prefix =~ s/\.gz.*//;
    $prefix =~ s/^\.\///;
    $prefix =~ s/\_.+//;


    my $alias="join";
    my $cmd_file=$out_dir."join.cmd";
    my $output_file = $out_dir.$prefix.".freqs";
    my $cmd1="perl $scripts_dir/join_freqs_file.pl $out_dir $output_file $do_gaps";

    Create_cmd::create_cmd_file($cmd_file,$alias,1,2,$cmd1);

    $cmd="qsub $cmd_file";

    my $stdin=`$cmd`;
    &sleep_manager(5,$stdin,$alias);
    unless (-e $output_file) {
	die "ERROR! Final output file was not created\n";
    }
}# join freq files





sub get_number_of_fasta_files_there_should_be {

    my $list_file=shift;   
    if ($num_fasta_files_in_this_run<0) {
	   $num_fasta_files_in_this_run=&get_num_fasta_files($list_file);

    }
    return $num_fasta_files_in_this_run;

}


sub get_num_fasta_files {

    print "get_num_fasta_files\n";

    my $list_file=shift;
    my $cmd="wc $list_file \| awk \'{print \$1}\'";
    my $num_files=`$cmd`;
    die "error ,get_num_fasta_files returned zero\n" if ($num_files==0);

    print "get num fasta END\n";
    return $num_files;
}



sub sleep_manager {

    my $stage_number=shift;
    my $stdin=shift;
    my $job_name=shift;

#4797045.power2.tau.ac.il

    die "unxpected error, cannot find job number, qsub returned $stdin\n" unless ($stdin =~ m/^(\d+).*/);
    my $job_num=$1;
    my $i=0;
    sleep(30);

    while ($i<$sleep_max) {
    	if (&test_qstat($job_num)>0)  {
    	    print "sleeping $i...";
    	    sleep ($i);
    	}
    	else {
    	    print "$job_name done, no more jobs in q\n";
    	    return;
    	}
	   $i+=$sleep_quantum;
    }

    sleep(10); # sleep another ten seconds to ensure that files appear for ls (mounting bug!! on qb3 cluster)

    print "\nexceeded max sleeping time $sleep_max\n";
}



sub test_qstat {

    my $job_num=shift;
    my $cmd = "qstat -t \| awk \'\$1 ~ \/^$job_num\/' \| wc \| awk \'{print \$1}\'";
    my $num_q_jobs=`$cmd`; chomp $num_q_jobs;
    
    return $num_q_jobs;
}





sub wrap_up {

    my $tmp_dir = $out_dir."tmp/";
    system("mkdir $tmp_dir") unless (-e $tmp_dir);

    my $cmd = "mv $out_dir"."*part* "."$tmp_dir";
    system($cmd);

    $cmd = "mv *.o* "."$tmp_dir";
    system($cmd);
   
    $cmd = "touch $tmp_dir/.ARK_NOBACKUP"; # indicated this dir is not backed up
    system($cmd);

    my $mutations_files = $tmp_dir.'*mutations.list';
    my $mutations_outfile = $out_dir.'mutations_all.txt'; 

    system("cat $mutations_files > $mutations_outfile");
}



