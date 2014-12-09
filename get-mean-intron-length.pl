#!/usr/bin/perl
use strict;

my $usage="$0 introns.fasta";
die "$usage\n" unless @ARGV==1;
my ($infile)=@ARGV;

my ($count,$sum);
open(IN,$infile) || die "can't open $infile\n";
while(<IN>)
  {
    if(/length=(\d+)/)
      {
	$sum+=$1;
	++$count;
      }
  }
close(IN);
my $mean=int(0.5+$sum/$count);
print "$mean\n";





