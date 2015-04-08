#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <chunks-dir>\n" unless @ARGV==1;
my ($dir)=@ARGV;

my @files=`ls $dir`;
my $n=@files;
for(my $i=0 ; $i<$n ; ++$i) {
  my $file=$files[$i];
  chomp $file;
  next unless $file=~/^(\d+)-alt.fasta/;
  my $id=$1;
  open(IN,"$dir/$file") || die "can't open file: $file";
  my $defline=<IN>; chomp $defline;
  close(IN);
  $defline=~/\/cigar=(\S+)/ || die "can't parse defline: $defline";
  my $cigar=$1;
  my $outfile="$dir/$id.lab";
  my $cmd="genezilla/project-annotation $dir/$id-ref.gff $dir/$id-ref.fasta $dir/$id-alt.fasta $cigar $outfile";
  system("$cmd");
  #print "$cmd\n";
}

