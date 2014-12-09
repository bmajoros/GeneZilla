#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;
use GffTranscriptReader;

my $maxLen=
my $separator="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
my $separatorLen=length $separator;
system("mkdir new-chunks");

my $ls=`ls chunks/*.gff`;
my @ls=split/\s+/,$ls;
my $bigSeq;
my $numFiles=@ls;
my $gffReader=new GffTranscriptReader;
open(OUT,">new-chunks/1.gff") || die;
for(my $i=0 ; $i<$numFiles ; ++$i)
  {
    my $file=$ls[$i];
    $file=~/(\d+)\.gff/ || die $file;
    my $filestem=$1;
    my $fastaIn="chunks/$filestem.fasta";
    my $gffIn="chunks/$filestem.gff";
    my $transcripts=$gffReader->loadGFF($gffIn);
    my $transcript=$transcripts->[0];
    my $fastaReader=new FastaReader($fastaIn);
    my ($def,$seq)=$fastaReader->nextSequence();
    $def=~/>(\S+)/ || die $def;
    my $seqId=$1;
    if($i>0) {$bigSeq.=$separator}
    $transcript->shiftCoords(length $bigSeq);
    $transcript->setSubstrate(1);
    my $gff=$transcript->toGff();
    print OUT "$gff\n";
    $bigSeq.=$seq;
    $fastaReader->close();
  }
close(OUT);
my $fastaWriter=new FastaWriter;
$fastaWriter->writeFasta(">1",$bigSeq,"new-chunks/1.fasta");






