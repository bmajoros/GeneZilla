#!/usr/bin/perl
use strict;
use ProgramName;

my $LOWEST=-1000;

my $name=ProgramName::get();
my $usage="$name <old.model> <new.model>";
die "$usage\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

open(IN,$infile) || die "Can't open $infile\n";
open(OUT,">$outfile") || die "Can't create $outfile\n";
while(<IN>)
  {
    if(/^([ATCGN]+)\s*$/)
      {
	print OUT;
	my $gram=$1;
	my $score=<IN>;
	if($gram=~/N/) {print OUT "$LOWEST\n"}
	else {print OUT "$score"}
      }
    else {print OUT}
  }
close(OUT);
close(IN);

