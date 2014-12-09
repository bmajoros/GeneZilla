#!/usr/bin/perl
use strict;

my $usage="$0 <approx-num-genes> <approx-gene-length> <genome-length>";
die "$usage\n" unless @ARGV==3;
my ($numGenes,$geneLength,$genomeSize)=@ARGV;

my $totalIntergenic=$genomeSize-$numGenes*$geneLength;
my $meanIntergenicSize=int($totalIntergenic/($numGenes+1)+0.5);

print "$meanIntergenicSize\n";



