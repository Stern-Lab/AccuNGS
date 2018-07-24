#!/usr/local/perl-5.20.2/bin/perl -w
use strict;
use lib '/usr/local/perl-5.20.2/lib/5.20.2/';

###################################################################################################

# Join numerous frequency files into one file, by multiplying freqs by number of reads.

###################################################################################################



die "usage join_freqs_files.pl <input files directory and prefix. suffix assumed is part-NN.fasta> <output file> <doGaps? Y/N default Y>" unless (scalar (@ARGV)>=2);


my @bases=("A","C","G","T","-");


my $input=$ARGV[0];
my $output_file = $ARGV[1];
my $do_gaps = "Y";
if (defined $ARGV[2]) {
    $do_gaps = $ARGV[2];
}
if ($do_gaps eq "N") {
    @bases=("A","C","G","T");
}

my %freqs;
my %sum_of_reads_for_pos;
my %ref_genome;

&main;

sub main {
    &analyse_files;
    &obtain_freqs_and_out;
}


sub obtain_freqs_and_out {
    open OUT, ">$output_file" or die "cannot open output file $output_file\n";
    print OUT "Pos\tBase\tFreq\tRef\tRead_count\tRank\tProb\n";

    foreach my $pos (sort {$a<=>$b} keys %freqs) {
	my @ranked_bases_order=&rank_and_test_freqs($pos);
	my $rank=0;
	foreach my $base_ind(@ranked_bases_order) {
	    my $base=$bases[$base_ind];
	    my $freq = sprintf "%.6f", $freqs{$pos}{$base};
	    my $prob = 1 - 10**((log(1-$freqs{$pos}{$base}+1e-07)/log(10))*($sum_of_reads_for_pos{$pos}{$base}+1));
	    $prob = sprintf "%.2f", $prob;
	    print OUT $pos."\t".$base."\t".$freq."\t".$ref_genome{$pos}."\t".$sum_of_reads_for_pos{$pos}{$base}."\t".$rank."\t".$prob."\n";
	    $rank++;
	}
    }
}

sub rank_and_test_freqs {
    my $pos=shift;
    my $sum=0;

    my @vals;
    foreach my $base (@bases) {
		
	if ($sum_of_reads_for_pos{$pos} ==0) {
	    $freqs{$pos}{$base}="?";
	}
	else {
	    $freqs{$pos}{$base}/=$sum_of_reads_for_pos{$pos}{$base};
	    push(@vals,$freqs{$pos}{$base});
	    $sum+=$freqs{$pos}{$base};
	    if ($freqs{$pos}{$base} <0 || $freqs{$pos}{$base} >1) {
		die "error, at position $pos, for base $base, obtained non-frequency value of ".$freqs{$pos}{$base}."\n";
	    }
	}
    }
    my $test=abs($sum-1.0);
    if ($test>0.01) {
	my $err_line=	 "ERROR, sum of freqs for pos $pos is $sum.\nTotal number of reads is".$sum_of_reads_for_pos{$pos}{"A"}. "\n";
	
	foreach my $base (@bases) {
	    $err_line.=$base.":".$freqs{$pos}{$base}."\n";
	}
	print $err_line."\n";
    }
    my @list_order = sort { $vals[$b] <=> $vals[$a] } 0 .. $#vals;

    return @list_order;   

}

sub analyse_files {
    my @files=glob($input."*part*.freqs");

    die "unexpected error, zero files match $input\n" if (scalar(@files)==0);
    foreach my $file(@files) {
	&analyse($file);
    }
}

sub analyse {
    my $file=shift;
    open IN, "<$file" or die "cannot open file $file\n";
    foreach my $line (<IN>) {
		next if ($line =~ m/\#/);
		my @f=split(/\t/,$line);
		my $pos=$f[0];
		my $base=$f[1];
		my $freq=$f[2];
		my $ref=$f[3];

		$ref_genome{$pos}=$ref;
		
		my $num_reads=$f[4];
		if ($freq eq "?") {
		    $freq=0;
		    $num_reads=0;
		}
		$sum_of_reads_for_pos{$pos}{$base}+=$num_reads;
		$freqs{$pos}{$base}+=$num_reads*$freq; # number of reads supporting this base call

    }
    close IN;

}
