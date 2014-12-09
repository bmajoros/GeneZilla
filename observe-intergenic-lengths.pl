#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use SummaryStats;

my $usage="$0 <*.gff>";
die "$usage\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $reader=new GffTranscriptReader();
my $transcripts=$reader->loadGFF($infile);

my $n=@$transcripts;
my @lengths;
my $prevTranscript;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    if(defined($prevTranscript) &&
       $prevTranscript->getSubstrate() eq $transcript->getSubstrate())
      {
	my $length=$transcript->getBegin()-$prevTranscript->getEnd();
	push @lengths,$length;
      }
    $prevTranscript=$transcript;
  }

my ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats(\@lengths);
print "mean=$mean, SD=$stddev, range=($min - $max)\n";

