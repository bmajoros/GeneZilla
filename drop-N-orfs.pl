#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;

my $usage="$0 <orfs.fasta> <out.fasta>";
die "$usage\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my $writer=new FastaWriter;
open(OUT,">$outfile") || die "can't create file: $outfile\n";
my $reader=new FastaReader($infile);
while(1)
  {
    my ($def,$seq)=$reader->nextSequence();
    last unless defined $def;
    next if $seq=~/N/;
    $writer->addToFasta($def,$seq,\*OUT);
  }
close(OUT);


