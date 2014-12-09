#!/usr/bin/perl
use strict;
use TempFilename;
use FastaReader;
use FastaWriter;

my $STANDARD_CONTEXT_LENGTH=80; # what the example-generator script provides

# PROCESS COMMAND LINE
$0=~/([^\/]+)$/;
my $usage=
"$1  <pos.fasta> <neg.fasta> <filestem> <\% train> 
     <signal-type> <consensus-offset> <consensus-length> 
     <context-window-length> <min-sens> <GC> <alpha> <min-sig-pos>

where
     <model-type> is one of {WMM,WAM,WWAM}
     <signal-type> is one of {ATG,TAG,GT,AG,PROMOTER,POLYA}
     <min-sig-pos> is the minimum number of significant positions
     <GC> is the genome's GC content, between 0 and 1 (ex: 0.48)
";
die "$usage\n" unless @ARGV==12;
my ($posFasta,$negFasta,$filestem,$percentTrain,$signalType,
    $consensusOffset,$consensusLength,$contextWindowLength,
    $minSensitivity,$GC,$alpha,$minSigPos)=@ARGV;

# DO SOME INITIALIZATION
my $genezilla=$ENV{"TIGRSCAN"};
if(!defined($genezilla)) 
  {die "environment variable TIGRSCAN must be set to GeneZilla install directory\n"}
my $dashG=($percentTrain<1 ? "-g" : "");

my $posExamples=substringFasta($posFasta,$consensusOffset,
			       $contextWindowLength,"pos-$signalType");
my $negExamples=substringFasta($negFasta,$consensusOffset,
			       $contextWindowLength,"neg-$signalType");

my $command="$genezilla/train-BWM $dashG ".
  "$posExamples $negExamples $filestem $percentTrain $signalType ".
  "$consensusOffset $consensusLength $minSensitivity $GC $alpha ".
  "$minSigPos";
print "$command\n";
system($command);


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





