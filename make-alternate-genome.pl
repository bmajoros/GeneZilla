#!/usr/bin/perl
use strict;
use FileHandle;
use FastaReader;
use FastaWriter;
use GffTranscriptReader;

die "\n
make-alternate-genome.pl <vcf-dir> <genome.2bit> <genes.gff> <out-dir>

Assumptions:
 * VCF files must be zipped with bgzip and have accompanying tbi indexes
 * VCF files must contain chromosome in filename, e.g. chr14, chrX, etc.
 * In VCF files, chrom names don't begin with \"chr\", but in GFF and
     2bit files, they do
 * twoBitToFa and tabix must be on your PATH
 * Your environment variable GENEZILLA must point to the GeneZilla dir
 * <out-dir> will be populated with two FASTA files per individual
"
  unless @ARGV==4;
my ($vcfDir,$twoBitFile,$gffFile,$outDir)=@ARGV;

#==============================================================
# First, some initialization
#==============================================================

my %chrToVCF;
initChrToVCF($vcfDir);
my $refGeneFasta="$outDir/refgene.fasta";
my $altGeneFasta="$outDir/altgene.fasta";
my $tempBedFile="$outDir/temp.bed";
my $geneVcfFile="$outDir/gene.vcf.gz";
my $geneGcfFile="$outDir/gene.gcf.gz";
my $GZ=$ENV{"GENEZILLA"};
my $fastaWriter=new FastaWriter;

#==============================================================
# Make FASTA files for each individual
#==============================================================

my $individuals=getIndividualList($vcfDir);
my $numIndiv=@$individuals;
my (@fastaFiles,%fastaFiles);
for(my $i=0 ; $i<$numIndiv ; ++$i) {
  my $indiv=$individuals->[$i];
  my $fh1=new FileHandle(">$outDir/$indiv-1.fasta");
  my $fh2=new FileHandle(">$outDir/$indiv-2.fasta");
  $fastaFiles{$indiv}=[$fh1,$fh2];
  push @fastaFiles,$fh1;
  push @fastaFiles,$fh2;
}
my $fh1=new FileHandle(">$outDir/ref-1.fasta");
my $fh2=new FileHandle(">$outDir/ref-2.fasta");
$fastaFiles{"reference"}=[$fh1,$fh2];
push @fastaFiles,$fh1; push @fastaFiles,$fh2;

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
  my $chr=$gene->getSubstrate();
  my $begin=$gene->getBegin();
  my $end=$gene->getEnd();
  my $strand=$gene->getStrand();
  my $name=$gene->getId();
  writeBed4($chr,$begin,$end,$name,$tempBedFile);
  System("twoBitToFa -bed=$tempBedFile -noMask $twoBitFile $refGeneFasta");
  my $chrVcfFile=$chrToVCF{$chr};
  writeBed3($chr,$begin,$end,$tempBedFile);
  System("tabix -h $chrVcfFile -R $tempBedFile | bgzip > $geneVcfFile");
  System("$GZ/vcf-to-gcf -c -v $geneVcfFile $geneGcfFile");
  writeBed6($chr,$begin,$end,$name,$strand,$tempBedFile);
  System("$GZ/gcf-to-fasta -r $geneGcfFile $twoBitFile $tempBedFile $altGeneFasta");
  my $fastaReader=new FastaReader($altGeneFasta);
  while(1) {
    my ($def,$seqref)=$fastaReader->nextSequenceRef();
    last unless $def;
    $def=~/>\S+\s+\/individual=(\S+)\s+\/allele=(\d+)\s+\/region=(\S+)/
      || die "Can't parse defline: $def\n";
    my ($indivID,$alleleNum,$geneID)=($1,$2,$3);
    my $fh=$fastaFiles{$indivID}->[$alleleNum];
    $def=">$geneID";
    $fastaWriter->addToFasta($def,$$seqref,$fh);
  }
  $fastaReader->close();
}

#==============================================================
# Clean up
#==============================================================

foreach my $file (@fastaFiles) { $file->close() }
unlink($refGeneFasta);
unlink($altGeneFasta);
unlink($tempBedFile);
unlink($geneVcfFile);
unlink($geneGcfFile);


#==============================================================
sub getIndividualList {
  my ($vcfDir)=@_;
  my @files=`ls $vcfDir/*.vcf.gz`;
  die "no VCF files found\n" unless @files>0;
  my $file=$files[0];
  my $individuals=[];
  open(IN,$file) || die "can't open file $file\n";
  while(<IN>) {
    chomp;
    if(/^\s*#CHROM/) {
      my @fields=split;
      my $numFields=@fields;
      for(my $i=9 ; $i<$numFields ; ++$i)
	{ push @$individuals,$fields[$i] }
      last;
    }
  }
  close(IN);
  return $individuals;
}
#==============================================================
# writeBed4($chr,$begin,$end,$name,$tempBedFile);
sub writeBed4 {
  my ($chr,$begin,$end,$name,$outfile)=@_;
  open(OUT,">$outfile") || die "Can't write file $outfile\n";
  print OUT "$chr\t$begin\t$end\t$name\n";
  close(OUT);
}
#==============================================================
# writeBed3($chr,$begin,$end,$tempBedFile);
sub writeBed3 {
  my ($chr,$begin,$end,$outfile)=@_;
  if($chr=~/chr(.+)/) { $chr=$1 }
  open(OUT,">$outfile") || die "Can't write file $outfile\n";
  print OUT "$chr\t$begin\t$end\n";
  close(OUT);
}
#==============================================================
# writeBed6($chr,$begin,$end,$name,$strand,$tempBedFile);
sub writeBed6 {
  my ($chr,$begin,$end,$name,$strand,$outfile)=@_;
  open(OUT,">$outfile") || die "Can't write file $outfile\n";
  print OUT "$chr\t$begin\t$end\t$name\t0\t$strand\n";
  close(OUT);
}
#==============================================================
sub System {
  my ($cmd)=@_;
  print "$cmd\n";
  system($cmd);
}
#==============================================================
# initChrToVCF($vcfDir);
sub initChrToVCF {
  my ($dir)=@_;
  my @files=`ls $dir/*.vcf.gz`;
  foreach my $file (@files) {
    chomp $file;
    if($file=~/(chr[A-Za-z\d]+)/) {
      my $chr=$1;
      $chrToVCF{$chr}=$file;
    }
  }
}
#==============================================================
#==============================================================
#==============================================================
#==============================================================
#==============================================================
#==============================================================
#==============================================================

