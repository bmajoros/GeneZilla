#!/usr/bin/perl
use strict;
use TempFilename;
use FastaReader;
use FastaWriter;

my $STANDARD_CONTEXT_LENGTH=80; # what the example-generator script provides

# PROCESS COMMAND LINE
$0=~/([^\/]+)\s*$/;
my $usage=
"$1 <model-type> <pos.fasta> <neg.fasta> <filestem> <\% train> 
     <signal-type> <consensus-offset> <consensus-length> 
     <context-window-length> <min-sens> <order> 
     <min-sample-size> <max-depth> <WWAM-wnd-size> <use-pooled-model(0/1)>";
die "$usage\n" unless @ARGV==15;
my ($modelType,$posFasta,$negFasta,$filestem,$percentTrain,$signalType,
    $consensusOffset,$consensusLength,$contextWindowLength,
    $minSensitivity,$order,$minSampleSize,$maxDepth,$WWAMwindowSize,
    $usePooledModel)=@ARGV;

# DO SOME INITIALIZATION
my $genezilla=$ENV{"GENEZILLA"};
if(!defined($genezilla)) 
  {die "environment variable GENEZILLA must be set to GeneZilla install directory\n"}
my $dashG=($percentTrain<1 ? "-g" : "");

my $posExamples=substringFasta($posFasta,$consensusOffset,
			       $contextWindowLength,"pos-$signalType");
my $negExamples=substringFasta($negFasta,$consensusOffset,
			       $contextWindowLength,"neg-$signalType");
my $dashP=($usePooledModel ? " -p " : "");
my $command="$genezilla/mdd $consensusOffset $consensusLength ".
  "$modelType $signalType $posExamples $negExamples $filestem ".
  "$minSensitivity $percentTrain -d $maxDepth -o $order -w $WWAMwindowSize ".
  "-s $minSampleSize $dashP";
print "$command\n";
system($command);
#unlink $posExamples;
#unlink $negExamples;


#---------------------------------------------------------------------
#---------------------------------------------------------------------
sub substringFasta
  {
    my ($infile,$consensusOffset,$contextWindowLength,$filestem)=@_;
    my $outfile="$filestem.tmp";#TempFilename::generate();
    print "reading \"$infile\"\n";
    my $reader=new FastaReader($infile);
    open(OUT,">$outfile") || die "can't create temp file $outfile";
    my $writer=new FastaWriter;
    my $begin=$STANDARD_CONTEXT_LENGTH-$consensusOffset;
    while(1)
      {
	my ($defline,$seq)=$reader->nextSequence();
	last unless defined $defline;
	my $subseq=substr($seq,$begin,$contextWindowLength);
	$writer->addToFasta($defline,$subseq,\*OUT);
      }
    close(OUT);
    return $outfile;
  }





