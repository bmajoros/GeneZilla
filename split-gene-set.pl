#!/usr/bin/perl
use strict;
use GffTranscriptReader;

die "split-gene-set.pl <in.gff> <num-subsets> <out-filestem>\n"
  unless @ARGV==3;
my ($infile,$numSets,$filestem)=@ARGV;

my $reader=new GffTranscriptReader;
my $genes=$reader->loadGenes($infile);
my $numGenes=@$genes;
my $genesPerSet=$numGenes/$numSets;
for(my $i=0 ; $i<$numSets ; ++$i) {
  my $outfile="$filestem$i.gff";
  my $begin=int($i*$genesPerSet);
  my $end=int(($i+1)*$genesPerSet);
  if($i==$numSets-1) { $end=$numGenes }
  open(OUT,">$outfile") || die "can't write to file $outfile";
  for(my $j=$begin ; $j<$end ; ++$j) {
    my $gene=$genes->[$j];
    my $gff=$gene->toGff();
    print OUT $gff;
  }
  close(OUT);
}

