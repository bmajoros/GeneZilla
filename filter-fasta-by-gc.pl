#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use FastaReader;
use FileHandle;
use FastaWriter;
use ProgramName;

my $name=ProgramName::get();
my $usage="$name <*.fasta> <min-GC> <max-GC>   [GC: decimal number from 0 to 1]";
die "$usage\n" unless @ARGV==3;
my ($infile,$minGC,$maxGC)=@ARGV;

my $fastaReader=new FastaReader($infile);
my $writer=new FastaWriter;
while(1)
  {
      my ($defline,$seqRef)=$fastaReader->nextSequenceRef();
      last unless defined $defline;
      my $GC=getGCcontent($seqRef);
      if($GC>=$minGC && $GC<$maxGC) 
      {
	  $writer->addToFasta($defline,$$seqRef,\*STDOUT);
      }
  }


#----------------------------------------------------------------
sub getGCcontent
  {
    my ($seqRef)=@_;
    my $GC=($$seqRef=~s/([GC])/$1/g);
    my $ATGC=($$seqRef=~s/([ATGC])/$1/g);
    $GC/=$ATGC;
    return $GC;
  }


