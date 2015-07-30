#!/usr/bin/perl
use strict;
use GffTranscriptReader;

my $MODEL="/home/bmajoros/splicing/model/no-UTR.iso";
my $NUC_MATRIX="/home/bmajoros/alignment/matrices/NUC.4.4";
my $GAP_OPEN=10;
my $GAP_EXTEND=1;
my $BANDWIDTH=100;

die "cia.pl <dir> <ID> <ref.fasta> <ref.gff> <alt.fasta>" unless @ARGV==5;
my ($dir,$ID,$refFasta,$refGff,$altFasta)=@ARGV;

#================================================================
# Do some initialization
#================================================================
my $GZ=$ENV{'GENEZILLA'};
if($GZ eq "") { die "Please set GENEZILLA environment variable\n" }
unless (`which fasta-seq-length.pl`=~/\/fasta-seq-length.pl/)
 { die "Can't find perl script fasta-seq-length.pl" }
my $cdsGff="$dir/$ID-cds.gff";
my $cigarFile="$dir/$ID.cigar";
my $labelFile="$dir/$ID.lab";
my $projectedGff="$dir/$ID-projected.gff";
my $variantsFile="$dir/$ID.variants";
my $signalsFile="$dir/$ID.signals";
my $graphFile="$dir/$ID.graph";
my $projectorReport="$dir/$ID-report.txt";

#================================================================
# First, check whether sequences have changed
#================================================================

my $changed=`$GZ/fastas-are-identical $refFasta $altFasta`;
if($changed eq "no") {
  print "No sequence changes -- terminating with no further analysis\n";
  exit(0);
}

#================================================================
# Exract CDS portion of gene model
#================================================================

System("grep CDS $refGff > $cdsGff");

#================================================================
# Reverse-complement files if necessary
#================================================================

my $refLen=0+`fasta-seq-length.pl $refFasta`;
my $strand=getStrand($cdsGff);
if($strand eq "-") {
  System("revcomp-gff.pl $cdsGff $refLen > $ID.tmp ; mv $ID.tmp $cdsGff");
  System("revcomp-fasta.pl $refFasta > $ID.tmp ; mv $ID.tmp $refFasta");
  System("revcomp-fasta.pl $altFasta > $ID.tmp ; mv $ID.tmp $altFasta");
}

#================================================================
# Align sequences to get CIGAR string
#================================================================

System("$GZ/BOOM/banded-smith-waterman -q -c $cigarFile $NUC_MATRIX $GAP_OPEN $GAP_EXTEND $refFasta $altFasta DNA $BANDWIDTH > /dev/null");

#================================================================
# Project annotation from ref to anno ("Blind Projector")
#================================================================

System("$GZ/project-annotation $cdsGff $refFasta $altFasta $cigarFile $labelFile $projectedGff");

#================================================================
# Check whether the projected gene model is broken
#================================================================

System("$GZ/check-projection $refFasta $cdsGff $altFasta $projectedGff $labelFile> $projectorReport");

#================================================================
# Find sequence variants and variant signals
#================================================================

System("$GZ/find-variants $refFasta $altFasta $cigarFile > $variantsFile");

System("$GZ/find-variant-signals $variantsFile $refFasta $cdsGff $altFasta $cigarFile $MODEL > $signalsFile");

#================================================================
# Build ORF graph using CIA model
#================================================================

System("$GZ/cia -O -g $graphFile $MODEL $altFasta $labelFile $projectedGff $signalsFile");

#================================================================
# Extract N best paths from graph
#================================================================

die;

System("$GZ/n-best");




#========================================================================
# SUBROUTINES
#========================================================================
sub System {
  my ($cmd)=@_;
  print "$cmd\n";
  system($cmd);
}
#========================================================================
sub getStrand {
  my ($gff)=@_;
  my $reader=new GffTranscriptReader();
  my $genes=$reader->loadGenes($gff);
  die "no genes in GFF file" unless @$genes>0;
  my $gene=$genes->[0];
  return $gene->getStrand();
}
#========================================================================
#========================================================================
#========================================================================
#========================================================================
#========================================================================
#========================================================================
#========================================================================
#========================================================================


