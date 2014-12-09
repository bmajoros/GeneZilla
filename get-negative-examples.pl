#!/usr/bin/perl
use strict;
$|=1;

my $usage="$0 <train-dir>";
die "$usage\n" unless @ARGV==1;
my ($dir)=@ARGV;

my $genezilla=$ENV{"GENEZILLA"};
if(!defined($genezilla)) 
  {die "environment variable GENEZILLA must be set to GeneZilla install directory\n"}

my $ls=`ls $dir/iso*-*.fasta`;
my @files=split/\s+/,$ls;
foreach my $fastaFile (@files)
  {
    $fastaFile=~/iso(.+)\.fasta/ || die;
    my $iso=$1;
    my $command="cd $dir ; $genezilla/get-negative-examples-iso.pl $iso";
    print "$command\nPlease wait, this may take a while...\n";
    system("$command");
  }

