#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use ProgramName;
my $name=ProgramName::get();
my $usage="$name <infile> <outfile>";
die "$usage\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my $reader=new GffTranscriptReader;
my $transcripts=$reader->loadGFF($infile);
open(OUT,">$outfile") || die "Can't write file: $outfile\n";
foreach my $transcript (@$transcripts)
  {
    my $substrate=$transcript->getSubstrate();
    if($substrate=~/^([^\s\.]+)\.(\S+)/)
      {
	my $geneId=$transcript->getGeneId();
	$transcript->setSubstrate($1);
	$transcript->setTranscriptId("$geneId.$2");
      }
    my $gff=$transcript->toGff();
    print OUT $gff;
  }
close(OUT);








