#!/usr/bin/perl
use strict;
use TempFilename;
use FastaReader;
use FastaWriter;

# PROCESS COMMAND LINE
my $usage="$0 <order> <min-sample-size> <model-type> <\%train> ".
"<content-type> <train-dir> <use-log-odds(0 or 1)>";
die "$usage\n" unless @ARGV==8;
my ($order,$minSampleSize,$modelType,$percentTrain,$contentType,
    $trainDir,$useLogOdds)=@ARGV;

# DO SOME INITIALIZATION
my $isochores=getIsochores($trainDir);
my $genezilla=$ENV{"TIGRSCAN"};
if(!defined($genezilla)) 
  {die "environment variable TIGRSCAN must be set to GeneZilla ".
     "install directory\n"}
my $filestem;
if($contentType eq "INTRON") {$filestem="introns"}
elsif($contentType eq "INTERGENIC") {$filestem="intergenic"}
elsif($contentType eq "SINGLE-EXON") {$filestem="single-exons"}
elsif($contentType eq "INITIAL-EXON") {$filestem="initial-exons"}
elsif($contentType eq "FINAL-EXON") {$filestem="final-exons"}
elsif($contentType eq "INTERNAL-EXON") {$filestem="internal-exons"}
else {die "That content type cannot be trained by this program.\n"}

# RUN TRAINER ON ALL ISOCHORES
my ($pooledPos,$pooledNeg);
foreach my $iso (@$isochores)
  {
    my $posExamples="$trainDir/$filestem$iso.fasta";
    my $negExamples="$trainDir/non-$filestem$iso.fasta";
    $pooledPos.="$posExamples ";
    $pooledNeg.="$negExamples ";
    my $nullModel="null$iso.model";
    my $dashR=($useLogOdds ? "-r$nullModel" : "");
    my $command="$genezilla/train-content-sensor $dashR -g -o $order -s $minSampleSize $modelType $posExamples $negExamples $filestem$iso $percentTrain $contentType";
    print "$command\n";
    system($command);
  }

# FOR COMPARISON, RUN TRAINING ON POOLED DATA
exit;
my $pooledPosFile=TempFilename::generate();
my $pooledNegFile=TempFilename::generate();
system("cat $pooledPos > $pooledPosFile");
system("cat $pooledNeg > $pooledNegFile");
my $command="$genezilla/train-content-sensor -g -o $order -s $minSampleSize $modelType $pooledPosFile $pooledNegFile $filestem-pooled $percentTrain $contentType";
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



