#!/usr/bin/perl
use strict;

my $usage="$0 <infile.trans> <outfile.trans>";
die "$usage\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

open(IN,$infile) || die "can't open $infile\n";
my (@lines,$total,$count);
open(OUT,">$outfile") || die "can't write into file: $outfile\n";
while(<IN>)
  {
    if(/(.*):\s*([\d\.\e]+)/)
      {
	push @lines,[$1,$2];
	$total+=$2;
	++$count;
      }
    else {print OUT}
  }
close(IN);
my $average=$total/$count;

foreach my $line (@lines)
  {
    my ($trans,$p)=@$line;
    $p/=$average;
    print OUT "$trans : $p\n";
  }
close(OUT);
