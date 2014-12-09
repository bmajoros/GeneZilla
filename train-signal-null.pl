#!/usr/bin/perl
use strict;
use TempFilename;
use FastaReader;
use FastaWriter;

my $STANDARD_CONTEXT_LENGTH=80; # what the example-generator script provides

# PROCESS COMMAND LINE
my $usage=
"$0 <model-type> <pos.fasta> <neg.fasta> <filestem> <\% train> 
     <signal-type> <consensus-offset> <consensus-length> 
     <context-window-length> <min-sens> <null-model> <order> 
     <min-sample-size>";
die "$usage\n" unless @ARGV==13;
my ($modelType,$posFasta,$negFasta,$filestem,$percentTrain,$signalType,
    $consensusOffset,$consensusLength,$contextWindowLength,
    $minSensitivity,$nullModel,$order,$minSampleSize)=@ARGV;

# DO SOME INITIALIZATION
my $genezilla=$ENV{"TIGRSCAN"};
if(!defined($genezilla)) 
  {die "environment variable TIGRSCAN must be set to GeneZilla install directory\n"}
my $dashG=($percentTrain<1 ? "-g" : "");
my $dashR="-r$nullModel ";

my $posExamples=substringFasta($posFasta,$consensusOffset,
			       $contextWindowLength,"pos-$signalType");
my $negExamples=substringFasta($negFasta,$consensusOffset,
			       $contextWindowLength,"neg-$signalType");
my $command="$genezilla/train-signal-sensor $dashR $dashG -o $order -s ".
  "$minSampleSize $modelType $posExamples ".
  "$negExamples $filestem $percentTrain $signalType ".
  "$consensusOffset $consensusLength $minSensitivity";
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





