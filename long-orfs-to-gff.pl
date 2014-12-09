#!/usr/bin/perl
use strict;

my $usage="
long-orfs-to-gff.pl <orfs> <substrate-id> <outfile>
   where <orfs> is the output of glimmer's long-orfs program";

die "$usage\n" unless @ARGV==3;
my ($infile,$contigId,$outfile)=@ARGV;

my $transcriptId=1;
open(IN,$infile) || die "Can't open $infile\n";
open(OUT,">$outfile") || die "Can't create $outfile\n";
while(<IN>)
  {
    if(/(\S+)\s+(\d+)\s+(\d+)/)
      {
	my ($id,$begin,$end)=($1,$2,$3);
	my $strand="+";
	if($begin>$end)
	  {
	    $strand="-";
	    ($begin,$end)=($end,$begin);
	  }
	print OUT "$contigId\tlong-orfs\tsingle-exon\t$begin\t$end\t.\t$strand\t0\ttransgrp=$transcriptId;\n";
	++$transcriptId;
      }
  }
close(OUT);
close(IN);



