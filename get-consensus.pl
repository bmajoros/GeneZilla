#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use FastaReader;

my $usage="$0 <isochore>     (e.g., 0-100)";
die "$usage\n" unless @ARGV==1;
my ($iso)=@ARGV;

my $MAX_EXAMPLES=2000;
my $DONOR_CONTEXT_BEGIN=-80;
my $DONOR_CONTEXT_LENGTH=162;
my $ACCEPTOR_CONTEXT_BEGIN=-80;
my $ACCEPTOR_CONTEXT_LENGTH=162;
my $START_CONTEXT_BEGIN=-80;
my $START_CONTEXT_LENGTH=163;
my $STOP_CONTEXT_BEGIN=-80;
my $STOP_CONTEXT_LENGTH=163;

getExamples("DONORS","donors$iso.fasta",80,2);
getExamples("ACCEPTORS","acceptors$iso.fasta",80,2);
getExamples("START-CODONS","start-codons$iso.fasta",80,3);
getExamples("STOP-CODONS","stop-codons$iso.fasta",80,3);

sub getExamples
  {
    my ($type,$filename,$contextLength,$consensusLength)=@_;
    my %examples;
    my $reader=new FastaReader($filename);
    print "\n$type\n";
    while(1)
      {
	my ($defline,$sequence)=$reader->nextSequence();
	last unless defined $defline;
	my $example=substr($sequence,$contextLength,$consensusLength);
	++$examples{$example};
      }

    my @examples=keys %examples;
    foreach my $consensus (@examples)
      {
	my $count=$examples{$consensus};
	print "$consensus\t$count\n";
      }
    print "\n";
  }
