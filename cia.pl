#!/usr/bin/perl
use strict;

die "cia.pl <ref.fasta> <ref.gff> <alt.fasta>" unless @ARGV==3;
my ($refFasta,$refGff,$altFasta)=@ARGV;

System("gz/fastas-are-identical gi_91176336_ref_NM_138383.2_.broken.fa  gi_91176336_ref_NM_138383.2_.fixed.fa");
#no

System("grep CDS gi_91176336_ref_NM_138383.2_.gff > CDS.gff");

System("fasta-seq-length.pl gi_91176336_ref_NM_138383.2_.fixed.fa");
#26847 bp

System("revcomp-gff.pl CDS.gff 26847 > tmp.gff ; mv tmp.gff CDS.gff");

System("revcomp-fasta.pl gi_91176336_ref_NM_138383.2_.fixed.fa > gi_91176336_ref_NM_138383.2_.fixed.fasta");

System("revcomp-fasta.pl gi_91176336_ref_NM_138383.2_.broken.fa > gi_91176336_ref_NM_138383.2_.broken.fasta");

System("gz/BOOM/banded-smith-waterman -q -c cigar/gi_91176336_ref_NM_138383.2_.cigar /home/bmajoros/alignment/matrices/NUC.4.4 10 1 gi_91176336_ref_NM_138383.2_.fixed.fasta  gi_91176336_ref_NM_138383.2_.broken.fasta  DNA 100 > /dev/null");
#7422M2I2M1D19422M

System("gz/project-annotation CDS.gff gi_91176336_ref_NM_138383.2_.fixed.fasta  gi_91176336_ref_NM_138383.2_.broken.fasta  cigar/gi_91176336_ref_NM_138383.2_.cigar  gi_91176336_ref_NM_138383.2_.lab  gi_91176336_ref_NM_138383.2_.projected.gff");

System("gz/check-projection  gi_91176336_ref_NM_138383.2_.fixed.fasta  CDS.gff    gi_91176336_ref_NM_138383.2_.broken.fasta   gi_91176336_ref_NM_138383.2_.projected.gff   gi_91176336_ref_NM_138383.2_.lab");

System("gz/find-variants  gi_91176336_ref_NM_138383.2_.fixed.fasta  gi_91176336_ref_NM_138383.2_.broken.fasta  cigar/gi_91176336_ref_NM_138383.2_.cigar  >  gi_91176336_ref_NM_138383.2_.variants");

System("gz/find-variant-signals  gi_91176336_ref_NM_138383.2_.variants   gi_91176336_ref_NM_138383.2_.fixed.fasta   CDS.gff   gi_91176336_ref_NM_138383.2_.broken.fasta    cigar/gi_91176336_ref_NM_138383.2_.cigar   ~/splicing/model/no-UTR.iso  >  signals.gff");

System("gz/cia -O -g gi_91176336_ref_NM_138383.2_.graph   ~/splicing/model/no-UTR.iso   gi_91176336_ref_NM_138383.2_.broken.fasta  gi_91176336_ref_NM_138383.2_.lab    gi_91176336_ref_NM_138383.2_.projected.gff   signals.gff");

System("gz/n-best");


#========================================================================
sub System {
  my ($cmd)=@_;
  print "$cmd\n";
  #system($cmd);
}
#========================================================================
#========================================================================
#========================================================================
#========================================================================
#========================================================================
#========================================================================
#========================================================================
#========================================================================
#========================================================================


