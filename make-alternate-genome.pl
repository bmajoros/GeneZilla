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
my $altGeneFasta="$outDir/altgene.fast";
my $tempBedFile="$outDir/temp.bed";
my $geneVcfFile="$outDir/gene.vcf.gz";
my $geneGcfFile="$outDir/gene.gcf.gz";
my $GZ=$ENV{"GENEZILLA"};

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
  #BED6: chr begin end name score strand
  writeBed6($chr,$begin,$end,$name,$strand,$tempBedFile);
  System("$GZ/gcf-to-fasta -r $geneGcfFile $twoBitFile $tempBedFile $altGeneFasta");

exit;

#  Now we have two haplotypes for many individuals in gene-alt.fasta.
#  Need to append to the two fasta files for each individual (in @fastaFiles)

}

#==============================================================
# Clean up
#==============================================================

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

