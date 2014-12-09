#!/usr/bin/perl
use strict;



die "THIS PROGRAM IS NOT YET IN WORKING CONDITION.  USE get-duration-distr.pl INSTEAD.\n";



my $windowSize=3;
my $binSize=30;

# PROCESS COMMAND LINE
my $usage="$0 <train-dir> <smoothing-iterations>";
die "$usage\n" unless @ARGV==2;
my ($trainDir,$iterations)=@ARGV;

my $isochores=getIsochores($trainDir);
my $genezilla=$ENV{"TIGRSCAN"};
if(!defined($genezilla)) 
  {die "environment variable TIGRSCAN must be set to GeneZilla install directory\n"}
foreach my $iso (@$isochores)
  {
    process("initial-exons$iso");
    process("internal-exons$iso");
    process("final-exons$iso");
    process("single-exons$iso");
  }

#---------------------------------------------------------------------
sub process
  {
    my ($filestem)=@_;
    my $filename="$trainDir/$filestem.fasta";
    my $outfile="$trainDir/$filestem.distr";
    my $command="$genezilla/get-duration-distr.pl $filename $binSize ".
      "$windowSize $iterations > $outfile";
    print "$command\n";
    system($command);
  }
#---------------------------------------------------------------------
sub getIsochores
  {
    my ($dir)=@_;
    my $ls=`ls $dir/iso*-*.fasta`;
    my @files=split/\s+/,$ls;
    my @isos;
    foreach my $fastaFile (@files)
      {
	$fastaFile=~/iso([^\/]+)\.fasta/ || die;
	my $iso=$1;
	push @isos,$iso;
      }
    return \@isos;
  }
#---------------------------------------------------------------------
