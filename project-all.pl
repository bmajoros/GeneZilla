#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <chunks-dir>\n" unless @ARGV==1;
my ($dir)=@ARGV;

my @files=`ls chunks`;
my $n=@files;
for(my $i=0 ; $i<$n ; ++$i) {
  my $file=$files[$i];
  chomp $file;
  next unless $file=~/^(\d+)-alt.fasta/;
  my $id=$1;
  
}

