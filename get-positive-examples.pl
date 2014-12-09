#!/usr/bin/perl
use strict;
$|=1;

my $usage="$0 <train-dir>";
die "$usage\n" unless @ARGV==1;
my ($dir)=@ARGV;

my $genezilla=$ENV{"TIGRSCAN"};
if(!defined($genezilla)) 
  {die "environment variable TIGRSCAN must be set to GeneZilla install directory\n"}

my $ls=`ls $dir/iso*-*.fasta`;
my @files=split/\s+/,$ls;
foreach my $fastaFile (@files)
  {
    $fastaFile=~/(.+)\.fasta/ || die;
    my $gffFile="$1.gff";
    my $command="cd $dir ; $genezilla/get-examples.pl $gffFile $fastaFile";
    print "$command\nPlease wait, this may take a while...\n";
    system("$command");
  }



