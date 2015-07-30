#!/usr/bin/perl
use strict;

die "make-alternate-genome.pl <variants.vcf> <genome.2bit> <regions.bed> <out.fasta>"
  unless @ARGV==4;
my ($vcfFile,$twoBitFile,$regeionsFile,$outFile)=@ARGV;


