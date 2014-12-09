#!/usr/bin/perl
use strict;
use TempFilename;
use FastaReader;
use FastaWriter;

my $STANDARD_CONTEXT_LENGTH=80; # what the example-generator script provides

# PROCESS COMMAND LINE
$0=~/([^\/]+)\s*$/;
my $usage=
"$1 <pos.fasta> <neg.fasta> <filestem> <\% train> <consensus-offset>
     <branch-context-length> <acceptor-context-length> <min-sens> 
     <branch-order> <acceptor-order> <min-sample-size> 
     <boosting-iterations> <boosting-percentile>
where
     <\% train> should be between 0 and 1
     <consensus-offset> is from the beginning of the entire sensor window
     <branch-context-length> is the size of the branch point WWAM
     <acceptor-context-length> is the size of the acceptor WAM
     <branch-order> and <acceptor-order> are for the WWAM & WAM, rspctvly
     <boosting-percentile> should be less than 0.2 for best results";
die "$usage\n" unless @ARGV==13;
my ($posFasta,$negFasta,$filestem,$percentTrain,$consensusOffset,
    $branchContextLength,$acceptorContextLength,$minSensitivity,
    $branchOrder,$acceptorOrder,$minSampleSize,$boostingIterations,
    $boostPercentile)=@ARGV;

# DO SOME INITIALIZATION
my $genezilla=$ENV{"GENEZILLA"};
if(!defined($genezilla)) 
  {die "environment variable GENEZILLA must be set to GeneZilla install directory\n"}
my $dashG=($percentTrain<1 ? "-g" : "");
my $dashB=($boostingIterations>0 ? "-b $boostingIterations" : "");
my $contextWindowLength=$branchContextLength+$acceptorContextLength;

# Extract positive & negative examples of just the sensor window
my $posExamples=substringFasta($posFasta,$consensusOffset,
			       $contextWindowLength,"pos-BranchAcceptors");
my $negExamples=substringFasta($negFasta,$consensusOffset,
			       $contextWindowLength,"neg-BranchAcceptors");

# Perform training
my $command="$genezilla/train-branch-acceptor $dashG -s $minSampleSize ".
  "$posExamples $negExamples $filestem $percentTrain ".
  "$branchContextLength $consensusOffset $minSensitivity ".
  "$branchOrder $acceptorOrder $dashB -p $boostPercentile";
print "$command\n";
system($command);

# Clean up
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





