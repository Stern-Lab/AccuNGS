#!/usr/local/perl-5.20.2/bin/perl -w
use strict;
use FileHandle;

$| = 1;

#########################################################################
# Parses Blast (-outfmt 6, btop format) of reference seqeunce against the reads
# and performs base-calling, based on >2 matches (HSPs) from the same read
# that support a character
# Input: blast output file in -m6 btop format
# Output: frequencies of each base at each position
#########################################################################

die "usage base_call_and_freqs.pl <blast output file in -m6 btop format format> <fastq (or fastq.gz) with quality scores> <reference sequence used with blast, fasta format> <output file for frequencies> <infer gaps? Y/N default Y> <min number repeats> <quality cutoff; default=30, requires all repeats to have this q or higher>" unless (scalar(@ARGV)>=4);

my $current_version_path = ""; #path to folder with AccuNGS scripts

my $infile = $ARGV[0];
my $fastq_file = $ARGV[1];

my $type_file = "";
if (($fastq_file =~ /\.fastq/) || ($fastq_file =~ /\.qual/)) {
    $type_file = "f";
}
if ($fastq_file =~ /\.fastq\.gz/) {
    $type_file  = "z";
}
die "error, second file must be either fastq or fastq.gz type file. got $fastq_file\n" if ($type_file eq "");

my $ref_fasta=$ARGV[2];
my $out_file = $ARGV[3];
if (-e $out_file) {print "WARNING, output file exists, will be overwritten\n"; sleep(5);}
my $summary_file = $out_file.".stats";
open SUMMARY, ">$summary_file" or die "cannot open file $summary_file\n";

my $do_gaps="N";  
if (defined $ARGV[4]) {
    $do_gaps=$ARGV[4];
}

# default quality cutoff per base is 30
my $min_for_calling = $ARGV[5];

my $quality_cutoff=30;

if (defined $ARGV[6]) {
   $quality_cutoff=$ARGV[6];
}

my $ascii_q_file = "$current_version_path/ascii_table_processed.txt";


my %relevant_reads;
my %ref;
my $max_read_length = 0;

my %reads_processed; # sanity check, used for filtering

my %read2qual; # after filtering, reads and their quality lines

my %mapped_reads; # counter for reads mapping several times; main object in this script
my %base_calls;
my %base_counts;
my @bases=("A","C","G","T","-");
my %reads_edges; # start and end of read; used to see if 85% of read map
my %asciiScores; # table with relevant lower Q score symbols

if ($do_gaps eq "N") {
    @bases=("A","C","G","T");
}
&main;

close SUMMARY;


sub main {
    &count_num_times_read_matches; # this is a filter that ensures that we look at reads that map more than twice to ref
    &read_ref_seq;
    print "processing reads\n";
    &process_reads();
    print "filtering\n";
    &filter();
    &read_ascii_quality_table;
    print "reading quality scores from fastq\n";
    &read_quality;
    print "quality and summary\n";
    &summarize();
    &print_and_out();

}


sub count_num_times_read_matches {
    my $cmd = "awk \'{print \$1}\' $infile \| uniq -c";
    my @res=`$cmd`;
    my %num_repeats;
    foreach my $elem (@res) {
	die "unexpected error\n" unless ($elem=~ m/(\d+)\s+(\S+)/);
	my $num=$1;
	my $read_id=$2;
	if ($num >=$min_for_calling) {	# TODO: in cirseq was set to be > 1 and >= min_for_calling in NGS 
	    $relevant_reads{$read_id}=1;
	}
	$num_repeats{$num}++;
    }
    print SUMMARY "Number of repeats in reads:\n";
    foreach my $key (sort {$a<=>$b} keys %num_repeats) {
	print SUMMARY $key."\t".$num_repeats{$key}."\n";
    }
}

sub read_ref_seq {
    open REF, "<$ref_fasta" or die "cannot open ref seq $ref_fasta\n";
    my $counter=1;
    foreach my $line (<REF>) {
	chomp $line;
	next if ($line =~ m/\>/);
	my @l=split(//,$line);
	foreach my $letter(@l) {
	    die "unxpected letter in sequence of $ref_fasta, non-word like, position $counter = $letter\n" unless ($letter =~ m/^[a-zA-Z]$/);
	    $letter =~ tr/[a-z]/[A-Z]/;
	    $ref{$counter}=$letter;
	    $counter++;
	}
    }
}

sub process_reads {
    open BLAST, "<$infile" or die "cannot open blast file $infile\n";
    foreach my $line (<BLAST>) {
	chomp $line;
	my @f=split(/\t/,$line);
	die "error in blast format, should be 8 fields\n" unless (scalar(@f)==8);
	my $read_id=$f[0];
	next unless (defined $relevant_reads{$read_id}); # skip reads that match only once (only one repeat)

	my $ref_start=$f[1];
	my $ref_end=$f[2];
	my $read_start=$f[3];
	my $read_end=$f[4];
	my $read_strand=$f[5];
	my $length=$f[6];
	my $aln=$f[7];
	
	# read_length approximation
	if ($max_read_length < $read_end){$max_read_length = $read_end;} 
	elsif ($max_read_length < $read_start){$max_read_length = $read_start;} 
	
	if ($ref_end<$ref_start ) {
	    die "line $line, ref end $ref_end is smaller than start $ref_start\n";
	}
	&map_read_edges($read_id,$read_start,$read_end);
	
	my ($aln_length,$mutations) = &parse_aln($aln); # returns a reference to a hash, with all mutations found in the read segment versus the ref
	die "unepxcted error, aln length from blast output is $length, but obtained $aln_length for line $line\n" unless ($length==$aln_length);
	
	my $position_to_print=$ref_start;
	my $position=$ref_start;
	my $read_position=$read_start;
	my $gap_counter=0;
	
	for my $i(1..$aln_length) {
	    $position_to_print=$position;
	    my $ref_letter=$ref{$position};
	    if (defined $mutations->{$i}) {
		my $pair_hash=$mutations->{$i};
		my $letter_at_ref_from_aln=$pair_hash->[0];
		die "error, position $position i=$i letter at ref is $ref_letter and in the aln it is $letter_at_ref_from_aln, line is $line\n" unless ($letter_at_ref_from_aln eq $ref_letter || $letter_at_ref_from_aln eq "-");
		my $letter_at_read=$pair_hash->[1];
		if ($letter_at_ref_from_aln eq "-") {
		    $gap_counter++;
		    $position--; # do not advance the counter
		    $position_to_print.=".$gap_counter";
		} elsif ($letter_at_read eq "-"){
			if ($read_strand eq "plus"){
				$read_position--;	#quality of the deletion would be the quality of the adjacent base
			} elsif ($read_strand eq "minus"){
				$read_position++;	#quality of the deletion would be the quality of the adjacent base
			} 
		}
		else {
		    $gap_counter=0;
		}
		 ## do not basecall a mutation (diff from ref) mapped at the end of the read
		$mapped_reads{$position_to_print}{$read_id}{$letter_at_read}{$read_position}=1;
		$reads_processed{$read_id}{$read_position}{$position_to_print}{$letter_at_read} = 1 if ($letter_at_read ne "-"); # no need to count read positions where there is a gap
	    }
	    else { # exact match to reference sequence
		$mapped_reads{$position_to_print}{$read_id}{$ref_letter}{$read_position}=1;
		$reads_processed{$read_id}{$read_position}{$position_to_print}{$ref_letter} = 1 if ($ref_letter ne "-");
	    }
	    $position++;
	    if ($read_strand eq "plus") {
		$read_position++;
	    }
	    elsif ($read_strand eq "minus") {
		$read_position--;
	    }
	    else {
		die "error, read polarity is undefined for line $line\n";
	    }

	}
    }
}

sub filter {

# Searching for invalid cases where a section of a read maps to two different places on the ref genome. This causes invalid base-calling.
# These cases are removed from the tally
    my $rem_double_map=0;
    foreach my $read_id (keys %reads_processed){

## filteration for double mapping
	foreach my $read_position (keys %{$reads_processed{$read_id}}) {
	    if (keys %{$reads_processed{$read_id}{$read_position}}>1) { # read maps to several locations
	    
		foreach my $pos_at_ref (keys %{$reads_processed{$read_id}{$read_position}}) {
		    foreach my $base (keys %{$reads_processed{$read_id}{$read_position}{$pos_at_ref}}) {
			undef %{$mapped_reads{$pos_at_ref}{$read_id}{$base}};

			$rem_double_map++;
		    }  
		}
	    }
	}
    }
    print SUMMARY $rem_double_map. " mapped bases removed due to double mapping\n";
}

sub read_quality {
    print "finding quality for ".scalar (keys %reads_processed)." reads\n";
    my $list_reads = $infile.".good_reads.list";
    my $list_reads_q = $infile.".good_reads.list.quality";
 
    my $grep_name="zgrep";
    if ($type_file eq "f") {$grep_name = "grep";} 
    
    open OUT, ">$list_reads" or die "cannot open file $list_reads\n";
    foreach my $read (keys %reads_processed) {
	print OUT "@".$read."\n"; # removed space before \n due to grep problems
    }
    
    my $cmd_line= "";
    if ($min_for_calling == 1){
    	$cmd_line= "$grep_name -F -A1 -f $list_reads $fastq_file \| grep -v \'^--$\' >$list_reads_q";
    }
    else{
    	$cmd_line= "$grep_name -F -A3 -f $list_reads $fastq_file \| grep -v \'^--$\' >$list_reads_q";
    }
	
    my $res=`$cmd_line`;
    
	#sleeps 1 minute if file does not exist
    sleep 60 if (!(-e $list_reads_q));
    
	print "opening file $list_reads_q\n";
    open IN, "<$list_reads_q" or die "cannot open file $list_reads_q\n";
    my @lines=<IN>;
    die "unexp error, number of lines returned from search **$list_reads_q** of quality should be even and it is ".scalar (@lines)."\n" if ((scalar (@lines)) %2!=0 || scalar (@lines)==0);

	for (my $i=0; $i<scalar(@lines)-1; $i=$i+2) {    
	my $line=$lines[$i];
	chomp $line;
	if ($line =~ m/\@(\S+)/) { 
	    my $read=$1;
	    my $q_line = $lines[$i+1];
	    $read2qual{$read}=$q_line;
	}
	else {
	    die "unexp error, line $line did not match header for quality scores\n";
	}
    }
	close (IN);
}

sub summarize {
    foreach my $pos (sort {$a <=> $b} keys %mapped_reads) {
	foreach my $base (@bases) { ## this will define that we are going over the bases we are interesetd in (ACGT with or without a gap)
	    $base_calls{$pos}{$base}=0;
	}
	$base_counts{$pos}=0;
    }
	
    # generate additional file for writing contributing mutations
    my $list_reads_mutations = $infile.".good_reads.mutations.list";
    open (CRO, ">$list_reads_mutations") || die "cannot open file $list_reads_mutations for writing\n";
    print CRO "ref_pos\tread\tbase\tread_position\n";
    
    my %reads_contributing;
    my %num_repeats;
    my %tmphash = map {$_=>1} @bases; # used only for hash search on line 315
    
    foreach my $pos (sort {$a <=> $b} keys %mapped_reads) {
	foreach my $read (keys %{$mapped_reads{$pos}}) {
	    LINE: foreach my $base (keys %{$mapped_reads{$pos}{$read}} ) {
		next unless (defined $tmphash{$base});
	## here: need to obtain read locations, and grep out the Q score; note we have to add a space for the search for grep otherwise many matches come up (e.f., 1443 will come up for 144)
		my $num_repeats=scalar keys %{$mapped_reads{$pos}{$read}{$base}};
		if ($num_repeats >= $min_for_calling ) {

		    my $avg_q_score=0; 

		    foreach my $read_position (keys %{$mapped_reads{$pos}{$read}{$base}}) {
			die "unxep error, could not find quality for read $read\n" unless (defined $read2qual{$read});

			last if ($read_position>=$max_read_length);
			last if ($read_position<0);
			next LINE unless ($num_repeats >= $min_for_calling) ;
			
			my $q_val = substr $read2qual{$read},$read_position-1,1;	# read position is 1 baed index (starting from 1 and not zero)
			if (! defined $asciiScores{$q_val}){
				print STDERR "pos=$pos,read=$read,base=$base,read2qual{read}=$read2qual{$read},read_position=$read_position,qval=$q_val\n";
			}

			$avg_q_score+=$asciiScores{$q_val};

		    } #foreach read position
		    if ($num_repeats<1) {
			$avg_q_score=0 ;
		    }
		    else {
			$avg_q_score/=$num_repeats;
		    }

			if ($avg_q_score>=$quality_cutoff) {
					$reads_contributing{$read}=1;
					$num_repeats{keys % {$mapped_reads{$pos}{$read}{$base}}}++;
					$base_calls{$pos}{$base}++;
					$base_counts{$pos}++;

					#prints all contributing mutation information
					if (($base ne '-') && ($pos !~ /\./) && ($base ne $ref{$pos})) { #mutation that is not a gap
						my @readPositions = keys % {$mapped_reads{$pos}{$read}{$base}};
						my $readPosStr = join(':',@readPositions);
						print CRO "$pos\t$read\t$base\t$readPosStr\n";
					}
				
		    } #$avg_q_score>=$quality_cutoff

		}  #if num_repeats>= $min_for_calling 
	    } #base
	} #read
    } #pos
    
    close (CRO);
    
    print SUMMARY scalar(keys %reads_contributing)." reads contributing to frequency counts\n";
    foreach my $key (keys %num_repeats) {
	print SUMMARY "$key repeats, num bases called=".$num_repeats{$key}."\n";
    }

} #summarize


sub print_and_out {
    open OUT, ">$out_file" or die "cannot open file $out_file\n";
    print OUT "#Pos\tBase\tFreq\tRef\tCoverage\tRank\n";
    foreach my $pos (sort {$a <=> $b} keys %base_calls) {
	if ($pos =~ m/\./) { # this represents an insertion in the reference seqeunce
	    next if ($base_counts{$pos}==0);
	    my $gap_freq=$base_calls{$pos}{"-"}/$base_counts{$pos};
	    next if ($gap_freq ==1);
	}

	foreach my $base (@bases) {
	    my $freq="?";
	    if ($base_counts{$pos} !=0) {
		$freq=$base_calls{$pos}{$base}/$base_counts{$pos};

	    }
	}

	my @ranked_bases_order=&rank_freqs($pos);
	my $rank=0;
	foreach my $base_ind(@ranked_bases_order) {
	    my $base=$bases[$base_ind];
	    my $freq=$base_calls{$pos}{$base}/$base_counts{$pos};
	    $freq = sprintf "%.6f", $base_calls{$pos}{$base};
	    my $ref_let="-";
	    if (defined $ref{$pos}) {$ref_let=$ref{$pos};}
	    print OUT $pos."\t".$base."\t".$freq."\t".$ref_let."\t".$base_counts{$pos}."\t".$rank."\n";
	    $rank++;
	}
    }
    close OUT;
}


sub rank_freqs {
    my $pos=shift;
    my $sum=0;
    
    my @vals;
    foreach my $base (@bases) {
	if ($base_counts{$pos} ==0) {
	    $base_calls{$pos}{$base}="?";
	}
	else {
	    $base_calls{$pos}{$base}/=$base_counts{$pos};
	    push(@vals,$base_calls{$pos}{$base});
	}
    }
    my @list_order = sort { $vals[$b] <=> $vals[$a] } 0 .. $#vals;

    return @list_order;   

}


sub parse_aln {
# for example: 149GNTN5TN26GN1GN1GNGNGN1TN8GA8GNAN2GNTNGN1ANTNAN18AN5AN3
    my $aln_str=shift;
    my $aln_length=shift;
    my @words=grep { /\S/ }  split(/\d+/,$aln_str);
    my @numbers=split(/\D+/,$aln_str);
    die "unexpected error, there should be one more number in the aln than the number of mutations, aln is $aln_str\n" if ((scalar(@words)+1)!=scalar(@numbers));
    my %mutations;
    my $pos=1; 
    for my $i(0..scalar(@words)-1) {
	$pos+=$numbers[$i];
	my $mut_str=$words[$i];
	die "unexpected error, aln $aln_str there are mutations that do not divide by two\n" if ((length $mut_str) %2 != 0);
	for (my $j=0; $j<length($mut_str); $j=$j+2) {
	    my $ref_letter=substr $mut_str,$j,1;
	    my $read_letter=substr $mut_str,$j+1,1;
	    $mutations{$pos}=[$ref_letter,$read_letter];
	    $pos++;
	}
    }
    my $report_aln_length=$pos+$numbers[-1]-1; # all the rest are matches to reference, this is the last identical position. Backing the coutner once since it advance at the end of the last for loop

    return $report_aln_length,\%mutations;
}


# maps the beginning and end of the the different mappings of the read
sub map_read_edges {
    my $id=shift;
    my $start=shift;
    my $end=shift;

    if ($start>$end) {
	my $temp=$start;
	$start=$end;
	$end=$temp;
    }
    $reads_edges{$id}{NUM_BASES}+=$end-$start+1;
    if (defined $reads_edges{$id}{START}) {
	$reads_edges{$id}{START}=$start if ($start<$reads_edges{$id}{START});
    }
    else {
	$reads_edges{$id}{START}=$start;
    }
    if (defined $reads_edges{$id}{END}) {
	$reads_edges{$id}{END}=$end if ($end>$reads_edges{$id}{END});
    }
    else {
	$reads_edges{$id}{END}=$end;
    }

}

sub read_ascii_quality_table {
    open (ASCII_HANDLE,$ascii_q_file);
    
    while (my $asciiLine = <ASCII_HANDLE>) {
	
	chomp $asciiLine;
	
	#0 33 !
	my @data = split(/\s+/,$asciiLine);
	my $score = $data[0];
	my $symbol = $data[2];
	$asciiScores{$symbol}=$score;

    }

    close (ASCII_HANDLE);

}

