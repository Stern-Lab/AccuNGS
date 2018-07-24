#!/usr/local/perl-5.20.2/bin/perl

package Create_cmd;
use strict;

sub create_cmd_file {
    my $n=scalar(@_);
    my $file_name=shift;
    my $alias=shift;
    my $num_jobs=shift;
    my $memory=shift;
    my @cmd_lines = @_;
    
    open OUT, ">$file_name" or die "create_cmd_file: cannot open file $file_name\n";
    print OUT "#!/bin/bash\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -r y\n";
    print OUT "#PBS -q ${queue_name}\n";

    print OUT "#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH \n";
    print OUT "#PBS -N $alias\n";
    print OUT "#PBS -l mem=$memory"."000mb\n";

if ($num_jobs>1) {
    print OUT "#PBS -J 1-$num_jobs\n\n";
}
print OUT "module load perl/perl-5.20.2\n";
    print OUT join("\n",@cmd_lines);
  close OUT;  
}

1;
