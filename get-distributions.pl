#!/usr/bin/perl
use strict;




die "OBSOLETE!  Use get-duration-distr.pl!\n";






# PROCESS COMMAND LINE
my $usage="$0 <train-dir> <smoothing-window-size> <smoothing-iterations>";
die "$usage\n" unless @ARGV==3;
my ($trainDir,$windowSize,$iterations)=@ARGV;

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
    my $command="$genezilla/get-distribution.pl $filename ".
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
