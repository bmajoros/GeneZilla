



OBSOLETE




#!/usr/bin/perl
use strict;
use TempFilename;
use FastaReader;
use FastaWriter;
use Getopt::Std;
our $opt_f;
getopts('f:');

if($opt_f) {print "f=$opt_f\n"}

my $STANDARD_CONTEXT_LENGTH=80; # what the example-generator script provides

# PROCESS COMMAND LINE
$0=~/([^\/]+)\s*$/;
my $usage="$1 <order> <min-sample-size> <model-type> <\%train> ".
"<signal-type> <consensus-offset> <consensus-length> ".
"<context-window-length> ".
"<min-sensitivity> <train-dir> <use-log-odds(0 or 1)>";
die "$usage\n" unless @ARGV==11;
my ($order,$minSampleSize,$modelType,$percentTrain,$signalType,
    $consensusOffset,$consensusLength,$contextWindowLength,
    $minSensitivity,$trainDir,$useLogOdds)
  =@ARGV;


# DO SOME INITIALIZATION
my $isochores=getIsochores($trainDir);
my $genezilla=$ENV{"TIGRSCAN"};
if(!defined($genezilla)) 
  {die "environment variable TIGRSCAN must be set to GeneZilla install directory\n"}
my $filestem;
if($signalType eq "ATG") {$filestem="start-codons"}
elsif($signalType eq "TAG") {$filestem="stop-codons"}
elsif($signalType eq "GT") {$filestem="donors"}
elsif($signalType eq "AG") {$filestem="acceptors"}
else {die "That signal type cannot be trained by this program.\n"}
my $dashG=($percentTrain<1 ? "-g" : "");
if($opt_f) {$filestem=$opt_f}

# RUN TRAINER ON ALL ISOCHORES
my ($pooledPos,$pooledNeg);
foreach my $iso (@$isochores)
  {
    my $posExamples=substringFasta("$trainDir/$filestem$iso.fasta",
				   $consensusOffset,$contextWindowLength,
				   "pos-$signalType$iso");
    my $negExamples=substringFasta("$trainDir/non-$filestem$iso.fasta",
				   $consensusOffset,$contextWindowLength,
				   "neg-$signalType$iso");
    $pooledPos.="$trainDir/$filestem$iso.fasta ";
    $pooledNeg.="$trainDir/non-$filestem$iso.fasta ";
    my $nullModel="null$iso.model";
    my $dashR=($useLogOdds ? "-r$nullModel" : "");
    my $command="$genezilla/train-signal-sensor $dashR $dashG -o $order -s ".
      "$minSampleSize $modelType $posExamples ".
	"$negExamples $filestem$iso $percentTrain $signalType ".
	  "$consensusOffset $consensusLength $minSensitivity";
    print "$command\n";
    system($command);
    #unlink $posExamples;
    #unlink $negExamples;
  }

# FOR COMPARISON, RUN TRAINING ON POOLED DATA
exit;
my $pooledPosFile=TempFilename::generate();
my $pooledNegFile=TempFilename::generate();
system("cat $pooledPos > $pooledPosFile");
system("cat $pooledNeg > $pooledNegFile");
my $posExamples=substringFasta($pooledPosFile,
			       $consensusOffset,$contextWindowLength,
			       "pos-$signalType-pooled");
my $negExamples=substringFasta($pooledNegFile,
			       $consensusOffset,$contextWindowLength,
			       "neg-$signalType-pooled");
my $command="$genezilla/train-signal-sensor $dashG -o $order -s ".
  "$minSampleSize $modelType $posExamples ".
  "$negExamples $filestem-pooled $percentTrain $signalType ".
  "$consensusOffset $consensusLength $minSensitivity";
print "$command\n";
system($command);
unlink $pooledPosFile;
unlink $pooledNegFile;

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





