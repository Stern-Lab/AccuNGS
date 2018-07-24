#!/usr/bin/perl

use strict;
use warnings;
use IO::Uncompress::Gunzip qw($GunzipError);

$|=1;


########Split a fastq file into N fasta files##########


die "usage <inFile .fastq.gz> <outFile prefix .fasta (will be divied into N files)> <X number of reads in N files, default X=50K)> <z for zipped files f for fastq files, default z>\n" unless (scalar(@ARGV)>=2);

my $fileName = $ARGV[0]; # file name with full path (e.g. /data/.../Sample_1/1_ATTACTCG_L001_R1_001.fastq.gz)
my $outFile = $ARGV[1]; # output file
my $num_reads_per_file = 70000;

if (defined $ARGV[2]) {
    $num_reads_per_file=$ARGV[2];
}

my $type_file="z"; # default - gzipped fastq files
if (defined $ARGV[3]) {
    $type_file=$ARGV[3];
    die "type_file must be either f or z" unless ($type_file eq "z" || $type_file eq "f");
}


&main;


sub main {
    my $res="";
    if ($type_file eq "z") {
		$res=`zcat <$fileName \| wc -l`;
		print "fasta split checkpoint 1.1";
    }
    else {

		$res = `cat $fileName \| wc -l`;
		print "fasta split checkpoint 1";
    }

    die "unexpected error, res line is empty\n" if ($res eq "");

    $res =~ s/^\s+|\s+$//g;
    print "***$res***\n";
    die "uenxpected error, number of lines ".$res." in file $fileName does not divide by 4\n" unless ($res % 4==0);

# divide into chunks

    my $num_reads=$res/4;

    print "number of reads for file is $num_reads, number of lines is ".$res."\n";

    my $num_files=int($num_reads/$num_reads_per_file)+1; # number of reads in chunk
    if (($num_reads % $num_reads_per_file) == 0){$num_files--;}

    print "number of files to split into =  $num_files, with $num_reads_per_file reads per file\n";


    my $chunk_size_lines=$num_reads_per_file*4; # number of lines in fastq files per chunk
    my $filehandle;

 	if ($type_file eq "z") {
		$filehandle = IO::Uncompress::Gunzip->new( $fileName ) or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
    }
    else { 
    	open($filehandle,$fileName) || die "cannot open $fileName\n";
    }

	my $next_text_input = "";
	my $rowNum = 0;
	my $file_counter = 1;
	my $curr_out_file = $outFile.".part$file_counter.fasta";

    unless ($file_counter>$num_files){
	    open (OUT,">",$curr_out_file) || die "cannot open $curr_out_file for writing\n";
	    # adding quality output file
		my $curr_qual_out_file = $outFile.".part$file_counter.qual";
	    open (OUTQ,">",$curr_qual_out_file) || die "cannot open $curr_qual_out_file for writing\n";
   }
    if ($type_file eq "z"){	    

	    while (! $filehandle->eof()){
		    for (1..4) {
            $next_text_input .= $filehandle->getline();
            $rowNum++;
        }

        
        my @split_text = split(/\n/, $next_text_input);
   		unless (scalar @split_text == 4) {
     	 	die "Error, fastQ entry doesn't have 4 lines: " . $next_text_input;
    	}

    	my ($name_line, $seq_line, $plus, $qual_line) = @split_text;
		unless ($name_line =~ /^\@/) { 
        	die "Error, cannot identify first line as read name line: " . $next_text_input;
    	}

    	print OUTQ $name_line."\n";
		print OUTQ $qual_line."\n";	
    	$name_line =~ s/^\@/>/;
    	print OUT $name_line."\n";
		print OUT $seq_line."\n";	

    	$next_text_input = "";

        if (($rowNum % $chunk_size_lines)==0){
			close (OUT); 
			close (OUTQ);

			$file_counter++;
			unless ($file_counter>$num_files){
				$curr_out_file = $outFile.".part$file_counter.fasta";
    			open (OUT,">",$curr_out_file) || die "cannot open $curr_out_file for writing\n";
    			my $curr_qual_out_file = $outFile.".part$file_counter.qual";
    			open (OUTQ,">",$curr_qual_out_file) || die "cannot open $curr_qual_out_file for writing\n";
			}
		}
	} #while

	    $filehandle->close()

    } # type_file eq "z"

    else {
	while (! eof($filehandle)) {		
        for (1..4) {
            $next_text_input .= <$filehandle>;
            $rowNum++;
        }
        my @split_text = split(/\n/, $next_text_input);
   		unless (scalar @split_text == 4) {
     	 	die "Error, fastQ entry doesn't have 4 lines: " . $next_text_input;
    	}   

    	my ($name_line, $seq_line, $plus, $qual_line) = @split_text;
		unless ($name_line =~ /^\@/) { 
        	die "Error, cannot identify first line as read name line: " . $next_text_input;
    	}

    

    	print OUTQ $name_line."\n";
		print OUTQ $qual_line."\n";	
    	
        $name_line =~ s/^\@/>/;

    	print OUT $name_line."\n";
		print OUT $seq_line."\n";

    	$next_text_input = "";

        if (($rowNum % $chunk_size_lines)==0){

			close (OUT);
			close(OUTQ);

			$file_counter++;
			$curr_out_file = $outFile.".part$file_counter.fasta";

    		open (OUT,">",$curr_out_file) || die "cannot open $curr_out_file for writing\n";
    		my $curr_qual_out_file = $outFile.".part$file_counter.qual";
    		open (OUTQ,">",$curr_qual_out_file) || die "cannot open $curr_qual_out_file for writing\n";
		}	
    	} #while
	} #else (type_file ne 'z')
} #main

