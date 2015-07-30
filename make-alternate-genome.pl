#!/usr/bin/perl
use strict;
use FileHandle;
use FastaReader;
use FastaWriter;
use GffTranscriptReader;

die "\n
make-alternate-genome.pl <vcf-dir> <genome.2bit> <genes.gff> <out-dir>

 * vcf files must be named *.vcf.gz
 * twoBitToFa must be on your PATH
 * your environment variable GENEZILLA must point to the GeneZilla dir
 * <out-dir> will be populated with two FASTA files per individual
"
  unless @ARGV==4;
my ($vcfDir,$twoBitFile,$gffFile,$outDir)=@ARGV;

#==============================================================
# First, some initialization
#==============================================================

my $refGeneFasta="$outDir/refgene.fasta";
my $tempBedFile="$outDir/temp.bed";

#==============================================================
# Make FASTA files for each individual
#==============================================================

my $individuals=getIndividualList($vcfDir);
my $numIndiv=@$individuals;
my @fastaFiles;
for(my $i=0 ; $i<$numIndiv ; ++$i) {
  my $indiv=$individuals->[$i];
  push @fastaFiles,new FileHandle(">$outDir/$indiv-1.fasta");
  push @fastaFiles,new FileHandle(">$outDir/$indiv-2.fasta");
}

#==============================================================
# Load gene coordinates from GFF file
#==============================================================

my $gffReader=new GffTranscriptReader();
my $genes=$gffReader->loadGenes($gffFile);

#==============================================================
# Process each gene
#==============================================================

my $numGenes=@$genes;
for(my $i=0 ; $i<$numGenes ; ++$i) {
  my $gene=$genes->[$i];
  my $begin=$gene->getBegin();
  my $end=$gene->getEnd();


"twoBitToFa -bed=gene.bed -noMask /data/reddylab/Reference_Data/hg19.2bit gene.fasta";


#  hg19 => twoBitToFa => gene-ref.fasta

#  tabix => gene.vcf => vcf-to-gcf => gene.gcf => gcf-to-fasta => gene-alt.fasta
#  (also takes gff file, outputs individual gene gff files)

#  Now we have two haplotypes for many individuals in gene-alt.fasta.
#  Need to concatenate into fasta files for single individuals: just open
#  two files per individual, make a single pass over gene-alt.fasta, write
#  each sequence into the appropriate file.

}

#==============================================================
# Clean up
#==============================================================

unlink($refGeneFasta);
unlink($tempBedFile);




#==============================================================
sub getIndividualList {
  my ($vcfDir)=@_;
  my @files=`ls $vcfDir/*.vcf.gz`;
  die "no VCF files found\n" unless @files>0;
  my $file=$files[0];
  open(IN,$file) || die "can't open file $file\n";
  while(<IN>) {
    if(/^\s*#CHROM/) {
      my @fields=split;
      my $individuals=[];
      my $numFields=@fields;
      for(my $i=9 ; $i<$numFields ; ++$i)
	{ push @$individuals,$fields[$i] }
    }
  }
  close(IN);
}
#==============================================================
#==============================================================
#==============================================================
#==============================================================
#==============================================================
#==============================================================
#==============================================================
#==============================================================
#==============================================================
#==============================================================

