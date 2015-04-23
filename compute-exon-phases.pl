#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.gff>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $reader=new GffTranscriptReader();
my $transcripts=$reader->loadGFF($infile);
my $totalLength=0;
foreach my $transcript (@$transcripts) {
  my $numExons=$transcript->numExons();
  for(my $i=0 ; $i<$numExons ; ++$i) {
    my $phase=$totalLength%3;
    print "$phase\t";
    my $exon=$transcript->getIthExon($i);
    my $exonLen=$exon->getLength();
    $totalLength+=$exonLen;
  }
  my $phase=$totalLength%3;
  print "final phase (should be zero) is: $phase\n";
}


