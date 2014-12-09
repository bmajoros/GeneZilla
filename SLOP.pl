#!/usr/bin/perl
################################################################
# SLOP.pl : Separate Local Optimization of Parameters
#           (only for use as a baseline for comparison)
#
# bmajoros@tigr.org
# 10/15/2004
################################################################
use strict;
use Getopt::Std;
$|=1;

my $genezillaDir=$ENV{"TIGRSCAN"};
die "Environment variable TIGRSCAN must be set to genezilla directory\n" 
  unless defined $genezillaDir;

our($opt_d);
getopts("d");

my $usage="$0 <logfile.txt>";
die "$usage\n" unless @ARGV==1;
my ($logFileName)=@ARGV;

#---------------------------------------------------------------
# Hard-coded parameters
#---------------------------------------------------------------
my $MIN_SAMPLE_SIZE=175;
my $SINGLE_EXON_POOLING_THRESHOLD=50; # we pool if fewer than this
my $MIN_NONCODING_MARGIN=10;  # min size for noncoding signal margin
my $MAX_NONCODING_MARGIN=45;  # max size for noncoding signal margin
my $MIN_CODING_MARGIN=0;      # min size for coding signal margin
my $MAX_CODING_MARGIN=6;      # max size for coding signal margin

#---------------------------------------------------------------
# Global variables
#---------------------------------------------------------------
my $config;             # hash table; contents of config file
my ($atg,$tag,$gt,$ag); # signal sensor training parms
my $baseline;           # struct containing baseline scores
my $dir=`pwd`;
$dir=~s/\s+//g;
my $numChunks;

#---------------------------------------------------------------
# Verify some things before we begin
#---------------------------------------------------------------
if(!-e "genezilla") 
  {die "please create a symbolic link to genezilla directory!\n"}
if(!-e "iso0-100.gff" || !-e "iso0-100.fasta")
  {die "iso0-100.fasta or iso0-100.gff is missing!\n"}
if(!-e 'chunks')
  {die "No chunks directory -- please use the make-test-set.pl and\n" .
     "combine-chunks scripts!\n"};
if(!-e 'out') {mkdir("out")}
if(!-e "initial-exons0-100.distr" || !-e "internal-exons0-100.distr" ||
   !-e "final-exons0-100.distr" || !-e "single-exons0-100.distr")
  {die "One or more *.distr files are missing: run get-duration-distr.pl\n"}

#---------------------------------------------------------------
# Miscellaneous initialization
#---------------------------------------------------------------
if(!-e "genezilla.top")
  {system("cp genezilla/polya0-100.model genezilla/TATA0-100.model " .
	  "genezilla/genezilla.top .")}
open(LOG,">$logFileName") || die "Can't creat file: $logFileName\n";
$numChunks=countChunks();
$config=makeDefaultConfig($dir);

#---------------------------------------------------------------
# Extract example sequences from training data, if necessary
#---------------------------------------------------------------
if(!-e "donors0-100.fasta")
  {
    print STDERR "Extracting examples of all features...\n";
    system("get-examples.pl iso0-100.gff iso0-100.fasta TAG,TGA,TAA " .
	   "notrim");
    system("get-negative-examples-iso.pl 0-100");
  }

#---------------------------------------------------------------
# Train the content sensors
#---------------------------------------------------------------
if(!-e "initial-exons0-100.model") 
  {print STDERR "Training content sensors...\n"}

if(!-e "initial-exons0-100.model")
  {system("genezilla/train-content-sensor -o 5 -s $MIN_SAMPLE_SIZE 3PIMM initial-exons0-100.fasta non-initial-exons0-100.fasta initial-exons0-100 1.0 INITIAL-EXON -g ; genezilla/compile-markov-chain initial-exons0-100.model")}

if(!-e "final-exons0-100.model") 
  {system("genezilla/train-content-sensor -o 5 -s $MIN_SAMPLE_SIZE 3PIMM final-exons0-100.fasta non-final-exons0-100.fasta final-exons0-100 1.0 FINAL-EXON -g ; genezilla/compile-markov-chain final-exons0-100.model")}

if(!-e "internal-exons0-100.model") 
  {system("genezilla/train-content-sensor -o 5 -s $MIN_SAMPLE_SIZE 3PIMM internal-exons0-100.fasta non-internal-exons0-100.fasta internal-exons0-100 1.0 INTERNAL-EXON -g ; genezilla/compile-markov-chain internal-exons0-100.model")}

if(!-e "single-exons0-100.model")
  {system("genezilla/train-content-sensor -o 5 -s $MIN_SAMPLE_SIZE 3PIMM single-exons0-100.fasta non-single-exons0-100.fasta single-exons0-100 1.0 SINGLE-EXON -g ; genezilla/compile-markov-chain single-exons0-100.model")}

if(!-e "introns0-100.model")
  {system("genezilla/train-content-sensor -o 5 -s $MIN_SAMPLE_SIZE IMM introns0-100.fasta non-introns0-100.fasta introns0-100 1.0 INTRON -g ; genezilla/compile-markov-chain introns0-100.model")}

if(!-e "intergenic0-100.model")
  {system("genezilla/train-content-sensor -o 5 -s $MIN_SAMPLE_SIZE IMM intergenic0-100.fasta non-intergenic0-100.fasta intergenic0-100 1.0 INTERGENIC -g ; genezilla/compile-markov-chain intergenic0-100.model")}

if(!-e "five-prime-UTR0-100.model")
  {system("genezilla/train-content-sensor -o 5 -s $MIN_SAMPLE_SIZE IMM five-prime-UTR0-100.fasta three-prime-UTR0-100.fasta five-prime-UTR0-100 1.0 FIVE-PRIME-UTR -g ; genezilla/compile-markov-chain five-prime-UTR0-100.model")}

if(!-e "three-prime-UTR0-100.model")
  {system("genezilla/train-content-sensor -o 5 -s $MIN_SAMPLE_SIZE IMM three-prime-UTR0-100.fasta five-prime-UTR0-100.fasta three-prime-UTR0-100 1.0 THREE-PRIME-UTR -g ; genezilla/compile-markov-chain three-prime-UTR0-100.model")}

if(!-e "single-exons0-100.model")
  {
    my $numSingleExons=`grep -c '>' single-exons0-100.fasta`;
    if($numSingleExons<$SINGLE_EXON_POOLING_THRESHOLD)
      {
	system("cat initial-exons0-100.fasta internal-exons0-100.fasta final-exons0-100.fasta single-exons0-100.fasta > pooled-exons0-100.fasta");
	system("genezilla/train-content-sensor -o 5 -s $MIN_SAMPLE_SIZE 3PIMM pooled-exons0-100.fasta non-internal-exons0-100.fasta single-exons0-100 1.0 SINGLE-EXON -g ; genezilla/compile-markov-chain single-exons0-100.model");
      }
    else
      {
	system("genezilla/train-content-sensor -o 5 -s $MIN_SAMPLE_SIZE 3PIMM single-exons0-100.fasta non-single-exons0-100.fasta single-exons0-100 1.0 SINGLE-EXON -g ; genezilla/compile-markov-chain single-exons0-100.model");
      }
  }

#---------------------------------------------------------------
# Get an initial transition distribution
#---------------------------------------------------------------
if(!-e "default.trans" || -z "default.trans")
  {
    print STDERR "Observing transition probabilities...\n";
    system("genezilla/get-trans-distr.pl iso0-100.gff > default.trans");
  }

#---------------------------------------------------------------
# Optimize the signal sensors independently
#---------------------------------------------------------------
getDefaultSignalSensorParms();
optimizeSignal($atg);
optimizeSignal($tag);
optimizeSignal($gt);
optimizeSignal($ag);

#---------------------------------------------------------------
# Optimize all remaining parameters independently as well
#---------------------------------------------------------------
$config->{"mean-intron-length"}=
  `get-mean-intron-length.pl introns0-100.fasta`;
$config->{"mean-intergenic-length"}=1000;
$config->{"mean-5'-UTR-length"}=50;
$config->{"mean-3'-UTR-length"}=50;

#---------------------------------------------------------------
# Perform a single evaluation (just to show how poor SLOP is)
#---------------------------------------------------------------
writeConfigFile($config,"slop");
my $scores=evaluate("slop");
printScore("final score ",$scores);
system("rm tmp*.*");
close(LOG);
system("evaluate-single-chunk.pl >> $logFileName");
# -------------------  END OF MAIN PROGRAM ---------------------




################################################################
#                          subroutines
################################################################
sub writeConfigFile
  {
    my ($config,$filestem)=@_;

    # Write the *.CFG file
    my $configFile="$dir/$filestem.cfg";
    open(OUT,">$configFile") 
      || die "can't write config file: $configFile\n";
    my @keys=keys %$config;
    foreach my $key (@keys)
      {
	my $value=$config->{$key};
	print OUT "$key\t= $value\n";
      }
    close(OUT);

    # Write the *.ISO file
    my $isoFile="$dir/$filestem.iso";
    open(OUT,">$isoFile") || die "Can't write file: $isoFile\n";
    print OUT "G+C <= 100% : $configFile\n";
    close(OUT);
  }
#===============================================================
sub getDefaultSignalSensorParms
  {
    $atg={};
    $atg->{type}="ATG";
    $atg->{filestem}="start-codons";
    $atg->{modelFilestem}="start-codons";
    $atg->{leftMargin}=20;
    $atg->{consensusLen}=3;
    $atg->{rightMargin}=3;
    $atg->{boostIterations}=0;
    $atg->{boostPercentile}=0.1;
    $atg->{wamOrder}=0;
    $atg->{sensitivity}=0.98;
    $atg->{configMember}="start-codon-model";

    $tag={};
    $tag->{type}="TAG";
    $tag->{filestem}="stop-codons";
    $tag->{modelFilestem}="stop-codons";
    $tag->{leftMargin}=3;
    $tag->{consensusLen}=3;
    $tag->{rightMargin}=20;
    $tag->{boostIterations}=0;
    $tag->{boostPercentile}=0.1;
    $tag->{wamOrder}=0;
    $tag->{sensitivity}=0.98;
    $tag->{configMember}="stop-codon-model";

    $gt={};
    $gt->{type}="GT";
    $gt->{filestem}="donors";
    $gt->{modelFilestem}="donors";
    $gt->{leftMargin}=3;
    $gt->{consensusLen}=2;
    $gt->{rightMargin}=20;
    $gt->{boostIterations}=0;
    $gt->{boostPercentile}=0.1;
    $gt->{wamOrder}=0;
    $gt->{sensitivity}=0.98;
    $gt->{configMember}="donor-model";

    $ag={};
    $ag->{type}="AG";
    $ag->{filestem}="acceptors";
    $ag->{modelFilestem}="acceptors";
    $ag->{leftMargin}=20;
    $ag->{consensusLen}=2;
    $ag->{rightMargin}=3;
    $ag->{boostIterations}=0;
    $ag->{boostPercentile}=0.1;
    $ag->{wamOrder}=0;
    $ag->{sensitivity}=0.98;
    $ag->{configMember}="acceptor-model";
  }
#===============================================================
sub trainSignal
  {
    my ($signalParms,$percentTrain,$refIteration)=@_;

    my $type=$signalParms->{type};
    my $filestem=$signalParms->{filestem};
    my $modelFilestem=$signalParms->{modelFilestem};
    my $leftMargin=$signalParms->{leftMargin};
    my $consensusLen=$signalParms->{consensusLen};
    my $rightMargin=$signalParms->{rightMargin};
    my $boostIterations=$signalParms->{boostIterations};
    my $boostPercentile=$signalParms->{boostPercentile};
    my $wamOrder=$signalParms->{wamOrder};
    my $sensitivity=$signalParms->{sensitivity};
    my $windowLen=$leftMargin+$consensusLen+$rightMargin;

    open(IN,"genezilla/train-signal.pl WAM ${filestem}0-100.fasta non-${filestem}0-100.fasta ${modelFilestem}0-100 $percentTrain $type $leftMargin $consensusLen $windowLen $sensitivity $wamOrder $MIN_SAMPLE_SIZE $boostIterations $boostPercentile |") || die "can't execute: train-signal.pl\n";
    my $accuracy;
    while(<IN>)
      {
	if(/Iteration \#(\S+):\s+accuracy=(\S+)/) 
	  {
	    my $acc=$2;
	    if($acc>$accuracy)
	      {
		$accuracy=$acc;
		if($refIteration) {$$refIteration=$1}
	      }
	  }
      }
    close(IN);
    return $accuracy;
  }
#===============================================================
# This function runs genezilla and evaluates the accuracy.  It
# assumes the caller has already written the *.cfg file.
sub evaluate
  {
    my ($filestem)=@_;

    # Run genezilla on test set
    my $isoFile="$filestem.iso";
    system("genezilla/run.pl $isoFile");

    # Run evaluate.pl to score the accuracy, echo to log file
    my $output;
    if($numChunks==1) {$output=`genezilla/evaluate-single-chunk.pl`}
    else {$baseline=`genezilla/evaluate.pl`}
    #print LOG $output;

    # Parse results
    $output=~/(\S+)\s*\|(\d+)%\s*(\d+)\s*%\s*(\d+)%\s*\|\s*(\d+)\%\s*(\d+)%\|(\d+)%\s*(\d+)%\s*\|\s*(\d+)%\s*(\d+)%\s*\|\s*(\d+)\s*%\s*\((\d+)\)\s*/ 
      || die "Can't parse evaluation output:\n$output\n";
    my ($program,$nucSens,$nucSpec,$nucAcc,$spliceSens,$spliceSpec,
	$stSens,$stSpec,$exonSens,$exonSpec,$exactGenes,$numGenes)=
	  ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12);

    # Return a struct containing scores
    my $scores={};
    $scores->{nucSens}=$nucSens;
    $scores->{nucSpec}=$nucSpec;
    $scores->{nucAcc}=$nucAcc;
    $scores->{spliceSens}=$spliceSens;
    $scores->{spliceSpec}=$spliceSpec;
    $scores->{stSens}=$stSens;
    $scores->{stSpec}=$stSpec;
    $scores->{exonSens}=$exonSens;
    $scores->{exonSpec}=$exonSpec;
    $scores->{geneAcc}=$exactGenes;
    $scores->{numGenes}=$numGenes;
    $scores->{exonAcc}=($exonSens+$exonSpec)/2; # technically wrong, but...
    return $scores;
  }
#===============================================================
sub printScore
  {
    my ($label,$scores)=@_;
    print LOG "$label [nuc=$scores->{nucAcc}% exon=$scores->{exonAcc}% " .
      "gene=$scores->{geneAcc}% ($scores->{numGenes} genes)]\n";
  }
#===============================================================
sub makeDefaultConfig
  {
    my ($dir)=@_;
    my $cfg={};
    $cfg->{"use-signal-thresholds"}="true";
    $cfg->{"signal-threshold-multiplier"}=1;
    $cfg->{"invert-signal-probabilities"}="false";
    $cfg->{"metamodel-topology"}="$dir/genezilla.top";
    $cfg->{"donor-model"}="$dir/donors0-100.model";
    $cfg->{"donor-consensus"}="GT";
    $cfg->{"acceptor-model"}="$dir/acceptors0-100.model";
    $cfg->{"acceptor-consensus"}="AG";
    $cfg->{"start-codon-model"}="$dir/start-codons0-100.model";
    $cfg->{"start-codon-consensus"}="ATG";
    $cfg->{"stop-codon-model"}="$dir/stop-codons0-100.model";
    $cfg->{"stop-codon-consensus"}="TAG|TGA|TAA";
    $cfg->{"polya-model"}="$dir/polya0-100.model";
    $cfg->{"polya-consensus"}="AATAAA|ATTAAA";
    $cfg->{"promoter-model"}="$dir/TATA0-100.model";
    $cfg->{"promoter-consensus"}="TATAAA|TATATA|CATAAA|CATATA";
    $cfg->{"initial-exons"}="$dir/initial-exons0-100.binmod";
    $cfg->{"internal-exons"}="$dir/internal-exons0-100.binmod";
    $cfg->{"final-exons"}="$dir/final-exons0-100.binmod";
    $cfg->{"single-exons"}="$dir/single-exons0-100.binmod";
    $cfg->{"initial-exon-lengths"}="$dir/initial-exons0-100.distr";
    $cfg->{"internal-exon-lengths"}="$dir/internal-exons0-100.distr";
    $cfg->{"final-exon-lengths"}="$dir/final-exons0-100.distr";
    $cfg->{"single-exon-lengths"}="$dir/single-exons0-100.distr";
    $cfg->{"mean-intron-length"}=0;
    $cfg->{"mean-intergenic-length"}=0;
    $cfg->{"mean-5'-UTR-length"}=50;
    $cfg->{"mean-3'-UTR-length"}=50;
    $cfg->{"introns"}="$dir/introns0-100.binmod";
    $cfg->{"intergenic"}="$dir/intergenic0-100.binmod";
    $cfg->{"3'-UTR"}="$dir/three-prime-UTR0-100.binmod";
    $cfg->{"5'-UTR"}="$dir/five-prime-UTR0-100.binmod";
    $cfg->{"probability-start-in-intron"}=0;
    $cfg->{"probability-start-in-intergenic"}=1;
    $cfg->{"probability-start-in-5pUTR"}=0;
    $cfg->{"probability-start-in-3pUTR"}=0;
    $cfg->{"transition-probabilities"}="$dir/default.trans";
    $cfg->{"queue-capacity"}=1;
    $cfg->{"exon-optimism"}=0;
    $cfg->{"intron-optimism"}=0;
    return $cfg;
  }
#===============================================================
sub countChunks
  {
    return 0+`ls chunks/*.fasta | wc -l`;
  }
#===============================================================
sub evalNoncodingMargin
  {
    my ($margin,$codingMargin,$bestMarginRef,$bestScoreRef,
	$noncodingMember,$parms)=@_;
    return if
      $margin<$MIN_NONCODING_MARGIN || $margin>$MAX_NONCODING_MARGIN;
    my $configMember=$parms->{configMember};
    my $modelFilestem="tmp$codingMargin-$margin";
    $parms->{modelFilestem}=$modelFilestem;
    $parms->{$noncodingMember}=$margin;
    my $score=trainSignal($parms,0.8);
    if($score>$$bestScoreRef)
      {
	$$bestMarginRef=$margin;
	$$bestScoreRef=$score;
	print LOG "\tMARGIN=$margin ACCUR=$score\n";
      }
  }
#===============================================================
sub evalCodingMargin
  {
    my ($margin,$noncodingMargin,$bestMarginRef,$bestScoreRef,
	$codingMember,$parms)=@_;
    return if $margin<$MIN_CODING_MARGIN || $margin>$MAX_CODING_MARGIN;
    my $configMember=$parms->{configMember};
    my $modelFilestem="tmp$margin-$noncodingMargin";
    $parms->{modelFilestem}=$modelFilestem;
    $parms->{$codingMember}=$margin;
    my $score=trainSignal($parms,0.8);
    if($score>$$bestScoreRef)
      {
	$$bestMarginRef=$margin;
	$$bestScoreRef=$score;
	print LOG "\tMARGIN=$margin ACCUR=$score\n";
      }
  }
#===============================================================
sub optimizeMargins
  {
    my ($parms,$codingMember,$noncodingMember)=@_;

    # Optimize the noncoding margin to within +/-5 bp
    my $configMember=$parms->{configMember};
    my $origParmFilestem=$parms->{filestem};
    my $codingMargin=$parms->{$codingMember};
    my ($bestMargin,$bestScore);
    for(my $margin=$MIN_NONCODING_MARGIN ; $margin<=$MAX_NONCODING_MARGIN ;
	$margin+=5)
      {	evalNoncodingMargin($margin,$codingMargin,\$bestMargin,\$bestScore,
			    $noncodingMember,$parms) }

    # Optimize the noncoding margin to within +/-1 bp
    for(my $margin=$bestMargin-4 ; $margin<$bestMargin ; ++$margin)
      {	evalNoncodingMargin($margin,$codingMargin,\$bestMargin,\$bestScore,
			    $noncodingMember,$parms) }
    for(my $margin=$bestMargin+1 ; $margin<$bestMargin+5 ; ++$margin)
      {	evalNoncodingMargin($margin,$codingMargin,\$bestMargin,\$bestScore,
			    $noncodingMember,$parms) }
    trainSignal($parms,1.0);

    $parms->{modelFilestem}=$origParmFilestem;
    $config->{$configMember}="$dir/${origParmFilestem}0-100.model";
    my $modelFilestem="tmp$codingMargin-$bestMargin";
    my $modelFilename="$dir/${modelFilestem}0-100.model";
    my $newname="$dir/${origParmFilestem}0-100.model";
    system("mv $modelFilename $newname");
    system("rm tmp*.model tmp*.scores");
    $parms->{$noncodingMember}=$bestMargin;

    # Optimize the coding margin by +/-1 increments
    my $noncodingMargin=$bestMargin;
    undef $bestMargin;
    undef $bestScore;
    for(my $margin=$MIN_CODING_MARGIN ; $margin<=$MAX_CODING_MARGIN ;
	++$margin)
      {	evalCodingMargin($margin,$noncodingMargin,\$bestMargin,\$bestScore,
			 $codingMember,$parms) }
    trainSignal($parms,1.0);

    $parms->{modelFilestem}=$origParmFilestem;
    $config->{$configMember}="$dir/${origParmFilestem}0-100.model";
    my $modelFilestem="tmp$bestMargin-$noncodingMargin";
    my $modelFilename="$dir/${modelFilestem}0-100.model";
    my $newname="$dir/${origParmFilestem}0-100.model";
    system("mv $modelFilename $newname");
    system("rm tmp*.model tmp*.scores");
    $parms->{$codingMember}=$bestMargin;
    return $bestScore;
  }
#===============================================================
sub optimizeSignal
  {
    my ($parms)=@_;
    my $signalType=$parms->{type};
    print LOG "Optimizing signal sensor: $signalType\n";
    my $filestem=$parms->{filestem};
    my $noncodingMargin=
      ($signalType=~/ATG|AG/ ? "leftMargin" : "rightMargin");
    my $codingMargin=
      ($noncodingMargin eq "leftMargin" ? "rightMargin" : "leftMargin");

    # Optimize margins
    my $baseScore=optimizeMargins($parms,$codingMargin,$noncodingMargin);

    # Optimize WAM order
    optimizeWamOrder($parms,\$baseScore);

    # Try boosting
    optimizeBoosting($parms,\$baseScore);
  }
#===============================================================
sub optimizeWamOrder
{
  my ($parms,$baseScoreRef)=@_;
  my $bestOrder=0;
  for(my $order=1 ; $order<=3 ; ++$order)
    {
      last unless tryWamOrder($order,$parms,$baseScoreRef);
      $bestOrder=$order;
    }
  trainSignal($parms,1.0);
}
#===============================================================
sub tryWamOrder
  {
    my ($order,$parms,$baseScoreRef)=@_;
    my $baseFilestem=$parms->{filestem};
    my $standardPath="$dir/${baseFilestem}0-100.model";
    my $configMember=$parms->{configMember};
    $parms->{modelFilestem}="tmp";
    my $originalOrder=$parms->{wamOrder};
    $parms->{wamOrder}=$order;
    my $score=trainSignal($parms,0.8);
    my $modelPath="$dir/tmp0-100.model";
    $config->{$configMember}=$standardPath;
    print LOG "\tTrying WAM order=$order...";
    if($score>$$baseScoreRef)
      {
	$$baseScoreRef=$score;
	system("mv $modelPath $standardPath");
	print LOG "success!  ACCUR=$score\n";
	return 1;
      }
    else {print LOG "no good, ACCUR drops to $score\n"}
    $parms->{wamOrder}=$originalOrder;
    system("rm $modelPath");
    return 0;
  }
#===============================================================
sub optimizeBoosting
{
  my ($parms,$baseScoreRef)=@_;
  my $bestIter=0;
  my $bestPercentile=0;
  for(my $ptile=0.60 ; $ptile>0 ; $ptile-=0.05)
    {
      my $iter=50;
      if(tryBoosting($iter,$ptile,$parms,$baseScoreRef))
	{
	  $bestPercentile=$ptile;
	  $bestIter=$parms->{boostIterations};
	}
    }
  $parms->{boostIterations}=$bestIter;
  $parms->{boostPercentile}=$bestPercentile;
  if($bestPercentile)
    {
      my $iterations=$parms->{boostIterations};
      my $ptile=int(100*$bestPercentile+0.5);
      print LOG "\t$iterations BOOST ITERATIONS \@ $ptile\%, ACCUR=$$baseScoreRef\n";
    }
  else
    {
      print LOG "\tSELECTED NO BOOSTING\n";
    }
  trainSignal($parms,1.0);
}
#===============================================================
sub tryBoosting
  {
    my ($numIter,$ptile,$parms,$baseScoreRef)=@_;
    my $baseFilestem=$parms->{filestem};
    my $standardPath="$dir/${baseFilestem}0-100.model";
    my $configMember=$parms->{configMember};
    $parms->{modelFilestem}="tmp";
    my $originalPtile=$parms->{boostPercentile};
    my $originalIterations=$parms->{boostIterations};
    $parms->{boostPercentile}=$ptile;
    $parms->{boostIterations}=$numIter;
    my $score=trainSignal($parms,1.0,\$numIter);
    my $modelPath="$dir/tmp0-100.model";
    $config->{$configMember}=$standardPath;
    if($score>$$baseScoreRef)
      {
	$parms->{boostIterations}=$numIter;
	trainSignal($parms,1.0);
	system("mv $modelPath $standardPath");
	$$baseScoreRef=$score;
	return $numIter;
      }
    $parms->{boostPercentile}=$originalPtile;
    $parms->{boostIterations}=$originalIterations;
    system("rm $modelPath");
    return 0;
  }
#===============================================================
