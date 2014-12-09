#!/usr/bin/perl
use strict;

for(my $i=1 ; $i<=100 ; ++$i)
{
    my $len=0+`fasta-seq-length.pl chunks/$i.fasta`;
    my $cmd="revcomp-gff.pl out/${i}rc.gff $len > out/${i}rcrc.gff";
    print "$cmd\n";
    system($cmd);
}

