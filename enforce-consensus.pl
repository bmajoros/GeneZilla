#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;
use TempFilename;

my $CONTEXT_BEGIN=-80;
my $CONTEXT_LENGTH=162;

my $usage="$0 <*.fasta> <consensus-list>   {overwrites input file!!}

example:
  $0 stop-codons.pl TAG,TGA,TAA
";
die "$usage\n" unless @ARGV==2;
my ($infile,$list)=@ARGV;

my @consensuses=split/,/,$list;
my $consensusLen=length $consensuses[0];
my %keep;
foreach my $consensus (@consensuses) {$keep{$consensus}=1}

my $reader=new FastaReader($infile);
my $writer=new FastaWriter;
my $tempFilename=TempFilename::generate();
open(OUT,">$tempFilename") || die "can't create $tempFilename\n";
while(1)
  {
    my ($def,$seq)=$reader->nextSequence();
    last unless defined $def;
    my $consensus=substr($seq,-$CONTEXT_BEGIN,$consensusLen);
    next unless $keep{$consensus};
    $writer->addToFasta($def,$seq,\*OUT);
  }
close(OUT);
system("mv $tempFilename $infile");





