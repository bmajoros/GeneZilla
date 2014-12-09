#!/usr/bin/perl
################################################################
# GRAPE : GRadient Ascent Parameter Estimation / for GeneZilla
# bmajoros@tigr.org
# 10/5/2004
#
# Results:
#   SLOP:  nuc=91% exon=62% gene=34% (162 genes)
#   GRAPE: nuc=96% exon=77.5% gene=57% (269 genes)
################################################################
use strict;
use Getopt::Std;
$|=1;

### NOTE TO SELF: this script can't handle extra splice sites...
###               nor MDD, nor isochores...

my $genezillaDir=$ENV{"GENEZILLA"};
die "Environment variable GENEZILLA must be set to genezilla directory\n"
  unless defined $genezillaDir;

our($opt_d,$opt_C,$opt_q,$opt_m);
getopts("dCq:m");
my $optimizeDurations=$opt_d;
my $queueCapacity=1;
if($opt_q>0) {$queueCapacity=$opt_q}


$0=~/([^\/]+)$/;
my $program=$1;
my $usage="
$program [options] <logfile.txt> <max-num-test-genes> <genome-GC-content> 
    where:
          -d = optimize exon duration distributions
               <genomic-GC-content> is between 0 and 1 (ex: 0.48)
          -C = don't combine chunks
          -q N = set noncoding queue capacity to N
          -m = train MLE models and quit

Prerequisites:
  1) Must have training set in ./iso0-100.gff (or other isochore boundaries)
  2) Must have test set in ./test.gff
  3) Must have contigs in ./iso0-100.fasta
  4) Run get-examples.pl iso0-100.gff iso0-100.fasta TAG,TGA,TAA notrim
  5) Run get-duration-distr.pl for each type of exon
     (and view results using xgraph to ensure smoothness)
";
die "$usage\n" unless @ARGV==3;
my ($logFileName,$maxTestSet,$GC)=@ARGV;

#---------------------------------------------------------------
# Hard-coded parameters
#---------------------------------------------------------------
my $MIN_SAMPLE_SIZE=175;
my $EXON_POOLING_THRESHOLD=50;# we pool if fewer than this
my $FORCE_EXON_POOLING=1;     # always pool exons?
my $MAX_BWM_SAMPLE_SIZE=40;   # use BWM if <40, otherwise use WAM
my $MIN_NONCODING_MARGIN=10;  # min size for noncoding signal margin
my $MAX_NONCODING_MARGIN=45;  # max size for noncoding signal margin
my $MIN_CODING_MARGIN=0;      # min size for coding signal margin
my $MAX_CODING_MARGIN=10;     # max size for coding signal margin

#---------------------------------------------------------------
# Global variables
#---------------------------------------------------------------
my $config;             # hash table; contents of config file
my ($atg,$tag,$gt,$ag); # signal sensor training parms
my $baseline;           # struct containing baseline scores
my $axisMap;            # maps axis indices to axis names & vice-versa
my $pointScoreCache={}; # for memoization of search space
my $transitions;        # transition probabilities
my $dir=`pwd`;
$dir=~s/\s+//g;
my $numChunks;
my $isochore=determineIsochore(); # "iso#-#"
my ($initialExonLengths,$internalExonLengths,$finalExonLengths,
    $singleExonLengths);

#---------------------------------------------------------------
# Verify some things before we begin
#---------------------------------------------------------------
if(!-e "iso$isochore.gff") {die "iso$isochore.gff is missing!\n"}
if(!-e "test.gff") {die "test.gff is missing!\n"}
if(!-e "iso$isochore.fasta") {die "iso$isochore.fasta is missing!\n"}
if(!-e 'out') {mkdir("out")}
if(!-e "initial-exons$isochore.distr" || !-e "internal-exons$isochore.distr" ||
   !-e "final-exons$isochore.distr" || !-e "single-exons$isochore.distr")
  {die "One or more *.distr files are missing: run get-duration-distr.pl\n"}

#---------------------------------------------------------------
# Miscellaneous initialization
#---------------------------------------------------------------
if(!-e "chunks")
  {
    system("genezilla/make-test-set.pl test.gff iso$isochore.fasta $maxTestSet 1000 TGA,TAA,TAG");
    if(!$opt_C)
      {
	system("genezilla/combine-chunks.pl");
	system("mv chunks chunks-all ; mv new-chunks chunks");
      }
  }
if(!-e "genezilla.top")
  {system("cp genezilla/polya0-100.model genezilla/TATA0-100.model " .
	  "genezilla/genezilla.top .")}
open(LOG,">$logFileName") || die "Can't creat file: $logFileName\n";
$numChunks=countChunks();
system("cp initial-exons$isochore.distr  stable-initial.distr");
system("cp internal-exons$isochore.distr stable-internal.distr");
system("cp final-exons$isochore.distr    stable-final.distr");
system("cp single-exons$isochore.distr   stable-single.distr");

#---------------------------------------------------------------
# Extract example sequences from training data, if necessary
#---------------------------------------------------------------
if(!-e "donors$isochore.fasta")
  {
    print STDERR "Extracting examples of all features...\n";
    system("genezilla/get-examples.pl iso$isochore.gff iso$isochore.fasta TAG,TGA,TAA " .
	   "notrim");
  }
if(!-e "non-donors$isochore.fasta")
  {
    system("genezilla/get-negative-examples-iso.pl $isochore");
  }

#---------------------------------------------------------------
# Train the content sensors
#---------------------------------------------------------------

if(!-e "pooled-exons$isochore.fasta")
  {
    system("cat initial-exons$isochore.fasta internal-exons$isochore.fasta " .
	   "final-exons$isochore.fasta single-exons$isochore.fasta > " .
	   "pooled-exons$isochore.fasta")
  }

if(!-e "exons$isochore.model")
  {
    system("genezilla/train-content-sensor -o 5 -s 1 " .
	   "3PIMM pooled-exons$isochore.fasta non-initial-exons$isochore.fasta " .
	   "exons$isochore 1.0 INITIAL-EXON -g ; " .
	   "genezilla/compile-markov-chain exons$isochore.model");
  }

if(!-e "introns$isochore.model")
  {system("genezilla/train-content-sensor -o 5 -s $MIN_SAMPLE_SIZE IMM introns$isochore.fasta non-introns$isochore.fasta introns$isochore 1.0 INTRON -g ; genezilla/compile-markov-chain introns$isochore.model")}

if(!-e "intergenic$isochore.model")
  {system("genezilla/train-content-sensor -o 5 -s $MIN_SAMPLE_SIZE IMM intergenic$isochore.fasta non-intergenic$isochore.fasta intergenic$isochore 1.0 INTERGENIC -g ; genezilla/compile-markov-chain intergenic$isochore.model")}

if(!-e "five-prime-UTR$isochore.model")
  {system("genezilla/train-content-sensor -o 5 -s $MIN_SAMPLE_SIZE IMM five-prime-UTR$isochore.fasta three-prime-UTR$isochore.fasta five-prime-UTR$isochore 1.0 FIVE-PRIME-UTR -g ; genezilla/compile-markov-chain five-prime-UTR$isochore.model")}

if(!-e "three-prime-UTR$isochore.model")
  {system("genezilla/train-content-sensor -o 5 -s $MIN_SAMPLE_SIZE IMM three-prime-UTR$isochore.fasta five-prime-UTR$isochore.fasta three-prime-UTR$isochore 1.0 THREE-PRIME-UTR -g ; genezilla/compile-markov-chain three-prime-UTR$isochore.model")}


#---------------------------------------------------------------
# Get an initial transition distribution
#---------------------------------------------------------------
if(!-e "default.trans")
  {
    print STDERR "Observing transition probabilities...\n";
    system("genezilla/get-trans-distr.pl iso$isochore.gff > default.trans");
  }
$transitions=loadTransProbs("default.trans");

#---------------------------------------------------------------
# Get a baseline evaluation
#---------------------------------------------------------------
print STDERR "Evaluating baseline parameterization...\n";
getDefaultSignalSensorParms();
trainSignal($atg);
trainSignal($tag);
trainSignal($gt);
trainSignal($ag);
$baseline=getBaseline();
printScore("Baseline score: ",$baseline);

exit if $opt_m;

#---------------------------------------------------------------
# Optimize the signal sensors independently
#---------------------------------------------------------------
optimizeSignal($atg);
optimizeSignal($tag);
optimizeSignal($gt);
optimizeSignal($ag);

#---------------------------------------------------------------
# Optimize all remaining parameters using gradient ascent
#---------------------------------------------------------------
$initialExonLengths=loadDistr($config->{"initial-exon-lengths"});
$internalExonLengths=loadDistr($config->{"internal-exon-lengths"});
$finalExonLengths=loadDistr($config->{"final-exon-lengths"});
$singleExonLengths=loadDistr($config->{"single-exon-lengths"});
my $initialExonStep=4/@$initialExonLengths;
my $internalExonStep=4/@$internalExonLengths;
my $finalExonStep=4/@$finalExonLengths;
my $singleExonStep=4/@$singleExonLengths;
gradientAscent();

#writeConfigFile($config,"grape");
if(-e "grape.cfg")
  {
    system("cp grape.cfg working.cfg"); ### I don't think this is necessary
  }
system("mv best.cfg grape.cfg");

open(OUT,">grape.iso") || die "Can't write file: grape.iso\n";
print OUT "G+C <= 100% : $dir/grape.cfg\n";
close(OUT);

system("rm tmp*.*");
close(LOG);
if($numChunks==1) {`genezilla/evaluate-single-chunk.pl >> $logFileName`}
else {`genezilla/evaluate.pl >> $logFileName`}
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
    my ($signalParms)=@_;

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
    my $numExamples=`grep -c '>' ${filestem}$isochore.fasta`;

    if($numExamples<$MAX_BWM_SAMPLE_SIZE)
      {
	# train a BWM to account for small sample size
	my $alpha=0.05;
	my $minSignificantPos=2;
	$sensitivity=1;
	system("genezilla/train-bwm.pl ${filestem}$isochore.fasta " .
	     "non-${filestem}$isochore.fasta ${modelFilestem}$isochore 1 " .
	       "$type $leftMargin $consensusLen $windowLen $sensitivity " .
	       "$GC $alpha $minSignificantPos");
	$signalParms->{BWM}=1;
      }
    else
      { # train a WAM instead
	system("genezilla/train-signal.pl WAM ${filestem}$isochore.fasta " .
	     "non-${filestem}$isochore.fasta ${modelFilestem}$isochore 1 " .
	       "$type $leftMargin $consensusLen $windowLen $sensitivity " .
	       "$wamOrder $MIN_SAMPLE_SIZE $boostIterations " .
	       "$boostPercentile");
      }
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
    else {$output=`genezilla/evaluate.pl`}
    #print LOG $output;

    # Parse results
    $output=~/(\S+)\s*\|\s*(\d+)\s*%\s*(\d+)\s*%\s*(\d+)\s*%\s*\|\s*(\d+)\s*%\s*(\d+)\s*%\s*\|\s*(\d+)\s*%\s*(\d+)\s*%\s*\|\s*(\d+)\s*%\s*(\d+)\s*%\s*(\d+)\s*%\s*\|\s*(\d+)\s*%\s*\((\d+)\)\s*/ 
      || die "Can't parse evaluation output:\n$output\n";
    my ($program,$nucSens,$nucSpec,$nucAcc,$spliceSens,$spliceSpec,
	$stSens,$stSpec,$exonSens,$exonSpec,$exonF,$exactGenes,$numGenes)=
	  ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13);

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
    $scores->{exonAcc}=$exonF;
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
# sub compareScores(A,B) returns:
#   -1 if A<B
#    0 if A=B
#   +1 if A>B
sub compareScores
  {
    my ($score1,$score2)=@_;

    # First, try to decide based on NUC ACC only
    my $nuc1=$score1->{nucAcc};
    my $nuc2=$score2->{nucAcc};
    if($nuc1<$nuc2) {return -1}
    if($nuc1>$nuc2) {return 1}

    # NUC ACC is equal, so look at exon & gene scores
    my $exonDiff=$score1->{exonAcc}-$score2->{exonAcc};
    my $geneDiff=$score1->{numGenes}-$score2->{numGenes};
    my $totalDiff=$exonDiff+$geneDiff;
    if($totalDiff<0) {return -1}
    if($totalDiff>0) {return 1}
    return 0;
  }
#===============================================================
sub getBaseline
  {
    # Get default config parms
    $config=makeDefaultConfig($dir);

    # Write config parms into *.cfg file
    my $filestem="baseline";
    writeConfigFile($config,$filestem);

    # Run genezilla and parse accuracy numbers
    return evaluate($filestem);
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
    $cfg->{"donor-model"}="$dir/donors$isochore.model";
    $cfg->{"donor-consensus"}="GT";
    $cfg->{"acceptor-model"}="$dir/acceptors$isochore.model";
    $cfg->{"acceptor-consensus"}="AG";
    $cfg->{"start-codon-model"}="$dir/start-codons$isochore.model";
    $cfg->{"start-codon-consensus"}="ATG";
    $cfg->{"stop-codon-model"}="$dir/stop-codons$isochore.model";
    $cfg->{"stop-codon-consensus"}="TAG|TGA|TAA";
    $cfg->{"polya-model"}="$dir/polya0-100.model";
    $cfg->{"polya-consensus"}="AATAAA|ATTAAA";
    $cfg->{"promoter-model"}="$dir/TATA0-100.model";
    $cfg->{"promoter-consensus"}="TATAAA|TATATA|CATAAA|CATATA";
    $cfg->{"initial-exons"}="$dir/exons$isochore.binmod";
    $cfg->{"internal-exons"}="$dir/exons$isochore.binmod";
    $cfg->{"final-exons"}="$dir/exons$isochore.binmod";
    $cfg->{"single-exons"}="$dir/exons$isochore.binmod";
    $cfg->{"initial-exon-lengths"}="$dir/initial-exons$isochore.distr";
    $cfg->{"internal-exon-lengths"}="$dir/internal-exons$isochore.distr";
    $cfg->{"final-exon-lengths"}="$dir/final-exons$isochore.distr";
    $cfg->{"single-exon-lengths"}="$dir/single-exons$isochore.distr";
    $cfg->{"mean-intron-length"}=100;###
    $cfg->{"mean-intergenic-length"}=1000;###
    $cfg->{"mean-5'-UTR-length"}=100;
    $cfg->{"mean-3'-UTR-length"}=100;
    $cfg->{"introns"}="$dir/intergenic$isochore.binmod";
    $cfg->{"intergenic"}="$dir/intergenic$isochore.binmod";
    $cfg->{"3'-UTR"}="$dir/intergenic$isochore.binmod";
    $cfg->{"5'-UTR"}="$dir/intergenic$isochore.binmod";
    $cfg->{"probability-start-in-intron"}=0;
    $cfg->{"probability-start-in-intergenic"}=1;
    $cfg->{"probability-start-in-5pUTR"}=0;
    $cfg->{"probability-start-in-3pUTR"}=0;
    $cfg->{"transition-probabilities"}="$dir/default.trans";
    $cfg->{"queue-capacity"}=$queueCapacity;
    $cfg->{"exon-optimism"}=1;
    $cfg->{"intron-optimism"}=0;
    $cfg->{"histogram-interpolation"}="true";
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
    trainSignal($parms);
    $config->{$configMember}="$dir/${modelFilestem}$isochore.model";
    writeConfigFile($config,"tmp");
    my $scoreStruct=evaluate("tmp");
    if(!defined($$bestScoreRef) || 
       compareScores($scoreStruct,$$bestScoreRef)>0)
      {
	$$bestMarginRef=$margin;
	$$bestScoreRef=$scoreStruct;
	printScore("\tnon-coding margin => $margin ",$scoreStruct);
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
    trainSignal($parms);
    $config->{$configMember}="$dir/${modelFilestem}$isochore.model";
    writeConfigFile($config,"tmp");
    my $scoreStruct=evaluate("tmp");
    if(!defined($$bestScoreRef) || 
       compareScores($scoreStruct,$$bestScoreRef)>0)
      {
	$$bestMarginRef=$margin;
	$$bestScoreRef=$scoreStruct;
	printScore("\tcoding margin => $margin ",$scoreStruct);
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

    $parms->{modelFilestem}=$origParmFilestem;
    $config->{$configMember}="$dir/${origParmFilestem}$isochore.model";
    my $modelFilestem="tmp$codingMargin-$bestMargin";
    my $modelFilename="$dir/${modelFilestem}$isochore.model";
    my $newname="$dir/${origParmFilestem}$isochore.model";
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

    $parms->{modelFilestem}=$origParmFilestem;
    $config->{$configMember}="$dir/${origParmFilestem}$isochore.model";
    my $modelFilestem="tmp$bestMargin-$noncodingMargin";
    my $modelFilename="$dir/${modelFilestem}$isochore.model";
    my $newname="$dir/${origParmFilestem}$isochore.model";
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

    if(!$parms->{BWM})
      {
	# Optimize WAM order
	optimizeWamOrder($parms,\$baseScore);
      }

    # Optimize sensitivity
    optimizeSensitivity($parms,\$baseScore);
	
    if(!$parms->{BWM})
      {
	# Try boosting
	optimizeBoosting($parms,\$baseScore);
      }
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
}
#===============================================================
sub tryWamOrder
  {
    my ($order,$parms,$baseScoreRef)=@_;
    my $baseFilestem=$parms->{filestem};
    my $standardPath="$dir/${baseFilestem}$isochore.model";
    my $configMember=$parms->{configMember};
    $parms->{modelFilestem}="tmp";
    my $originalOrder=$parms->{wamOrder};
    $parms->{wamOrder}=$order;
    trainSignal($parms);
    my $modelPath="$dir/tmp$isochore.model";
    $config->{$configMember}=$modelPath;
    writeConfigFile($config,"tmp");
    my $scoreStruct=evaluate("tmp");
    $config->{$configMember}=$standardPath;
    if(compareScores($scoreStruct,$$baseScoreRef)>0)
      {
	system("mv $modelPath $standardPath");
	$$baseScoreRef=$scoreStruct;
	printScore("\tWAM order => $order ",$scoreStruct);
	return 1;
      }
    $parms->{wamOrder}=$originalOrder;
    system("rm $modelPath");
    return 0;
  }
#===============================================================
sub optimizeSensitivity
{
  my ($parms,$baseScoreRef)=@_;
  my $bestSens=0.98;
  for(my $sens=0.99 ; $sens<=1 ; $sens+=0.01)
    {
      last unless trySensitivity($sens,$parms,$baseScoreRef);
      $bestSens=$sens;
    }
  print LOG "BEST SENS=$bestSens\n";###
}
#===============================================================
sub trySensitivity
  {
    my ($sens,$parms,$baseScoreRef)=@_;
    my $baseFilestem=$parms->{filestem};
    my $standardPath="$dir/${baseFilestem}$isochore.model";
    my $configMember=$parms->{configMember};
    $parms->{modelFilestem}="tmp";
    my $originalSens=$parms->{sensitivity};
    $parms->{sensitivity}=$sens;
    trainSignal($parms);
    my $modelPath="$dir/tmp$isochore.model";
    $config->{$configMember}=$modelPath;
    writeConfigFile($config,"tmp");
    my $scoreStruct=evaluate("tmp");
    $config->{$configMember}=$standardPath;
    printScore("SENS: $sens ",$scoreStruct);
    if(compareScores($scoreStruct,$$baseScoreRef)>0)
      {
	system("mv $modelPath $standardPath");
	$$baseScoreRef=$scoreStruct;
	printScore("\tsignal sensitivity => $sens ",$scoreStruct);
	return 1;
      }
    $parms->{sensitivity}=$originalSens;
    system("rm $modelPath");
    return 0;
  }
#===============================================================
sub optimizeBoosting
{
  my ($parms,$baseScoreRef)=@_;
  my ($bestPercentile,$bestIter);
  for(my $ptile=0.20 ; $ptile>0 ; $ptile-=0.10)
    {
      foreach my $iter (1,3,5,10,20)
	{
	  if(tryBoosting($iter,$ptile,$parms,$baseScoreRef))
	    {
	      $bestPercentile=$ptile;
	      $bestIter=$iter;
	    }
	}
    }
  if($bestPercentile)
    {
      my $iterations=$parms->{boostIterations};
      my $ptile=int(100*$bestPercentile+5/9);
      print LOG "SELECTED $iterations BOOSTING ITERATIONS AT $ptile\%\n";###
    }
  else
    {
      print LOG "SELECTED NO BOOSTING\n";
    }
}
#===============================================================
sub tryBoosting
  {
    my ($numIter,$ptile,$parms,$baseScoreRef)=@_;
    my $baseFilestem=$parms->{filestem};
    my $standardPath="$dir/${baseFilestem}$isochore.model";
    my $configMember=$parms->{configMember};
    $parms->{modelFilestem}="tmp";
    my $originalPtile=$parms->{boostPercentile};
    my $originalIterations=$parms->{boostIterations};
    $parms->{boostPercentile}=$ptile;
    $parms->{boostIterations}=$numIter;
    trainSignal($parms);
    my $modelPath="$dir/tmp$isochore.model";
    $config->{$configMember}=$modelPath;
    writeConfigFile($config,"tmp");
    my $scoreStruct=evaluate("tmp");
    $config->{$configMember}=$standardPath;
    printScore("BOOSTING: $numIter ITERATIONS AT $ptile ",
	       $scoreStruct);
    if(compareScores($scoreStruct,$$baseScoreRef)>0)
      {
	system("mv $modelPath $standardPath");
	$$baseScoreRef=$scoreStruct;
	printScore("\tboosting => $numIter @ $ptile ",$scoreStruct);
	return 1;
      }
    $parms->{boostPercentile}=$originalPtile;
    $parms->{boostIterations}=$originalIterations;
    system("rm $modelPath");
    return 0;
  }
#===============================================================
sub initAxisMap
  {
    my $hash={};

    # Add the axes in random order
    my @elems=("intronLength","intergenicLength","utrLength","optimism",
	      "ATG->GT","AG->GT","polyA","promoter");
    if($optimizeDurations)
      {
	push @elems,"initialExonLengths";
	push @elems,"internalExonLengths";
	push @elems,"finalExonLengths";
	push @elems,"singleExonLengths";
      }
    for(my $i=0 ; $i<@elems ; ++$i)
      {
	my $j=int(rand(@elems-$i-1))+$i;
	my $tmp=$elems[$i];
	$elems[$i]=$elems[$j];
	$elems[$j]=$tmp;
      }
    for(my $i=0 ; $i<@elems ; ++$i)
      {addAxisMapElem($elems[$i],$hash)}
    return $hash;
  }
#===============================================================
sub addAxisMapElem
  {
    my ($name,$hash)=@_;
    my @keys=keys %$hash;
    my $id=@keys/2;
    $hash->{$name}=$id;
    $hash->{$id}=$name;
  }
#===============================================================
sub setMember
  {
    my ($array,$name,$value)=@_;
    my $index=$axisMap->{$name};
    $array->[$index]=$value;
  }
#===============================================================
sub getMember
  {
    my ($array,$name)=@_;
    my $index=$axisMap->{$name};
    return $array->[$index];
  }
#===============================================================
sub initState
  {
    my $state=[];
    setMember($state,"intronLength",
	      {
	       stepSize=>100,
	       direction=>0,
	       minStepSize=>5,
	       min=>1,
	       max=>1000000,
	       type=>"int"
	      });
    setMember($state,"intergenicLength",
	      {
	       stepSize=>1000,
	       direction=>0,
	       minStepSize=>50,
	       min=>1,
	       max=>1000000,
	       type=>"int"
	      });
    setMember($state,"utrLength",
	      {
	       stepSize=>100,
	       direction=>0,
	       minStepSize=>5,
	       min=>1,
	       max=>1000000,
	       type=>"int"
	      });
    setMember($state,"optimism",
	      {
	       stepSize=>1,
	       direction=>0,
	       minStepSize=>0.05,
	       min=>-5,
	       max=>5,
	       type=>"float"
	      });

    setMember($state,"ATG->GT",
	      {
	       stepSize=>0.1,
	       direction=>0,
	       minStepSize=>0.003,
	       min=>0,
	       max=>1,
	       type=>"float"
	      });
    setMember($state,"AG->GT",
	      {
	       stepSize=>0.2,
	       direction=>0,
	       minStepSize=>0.003,
	       min=>0,
	       max=>1,
	       type=>"float"
	      });
    setMember($state,"polyA",
	      {
	       stepSize=>0.2,
	       direction=>0,
	       minStepSize=>0.003,
	       min=>0,
	       max=>1,
	       type=>"float"
	      });
    setMember($state,"promoter",
	      {
	       stepSize=>0.2,
	       direction=>0,
	       minStepSize=>0.003,
	       min=>0,
	       max=>1,
	       type=>"float"
	      });
    if($optimizeDurations)
      {
	setMember($state,"initialExonLengths",
		  {
		   stepSize=>[$initialExonStep,30],
		   direction=>[0,0],
		   minStepSize=>[$initialExonStep/10,1],
		   min=>[0,0],
		   max=>[1,0],
		   type=>"distribution",
		   tempName=>"tmp-initial.distr",
		   stableName=>"grape-initial.distr",
		   configMember=>"initial-exon-lengths"
		  });
	setMember($state,"internalExonLengths",
		  {
		   stepSize=>[$internalExonStep,30],
		   direction=>[0,0],
		   minStepSize=>[$internalExonStep/10,1],
		   min=>[0,0],
		   max=>[1,0],
		   type=>"distribution",
		   tempName=>"tmp-internal.distr",
		   stableName=>"grape-internal.distr",
		   configMember=>"internal-exon-lengths"
		  });
	setMember($state,"finalExonLengths",
		  {
		   stepSize=>[$finalExonStep,30],
		   direction=>[0,0],
		   minStepSize=>[$finalExonStep/10,1],
		   min=>[0,0],
		   max=>[1,0],
		   type=>"distribution",
		   tempName=>"tmp-final.distr",
		   stableName=>"grape-final.distr",
		   configMember=>"final-exon-lengths"
		  });
	setMember($state,"singleExonLengths",
		  {
		   stepSize=>[$singleExonStep,30],
		   direction=>[0,0],
		   minStepSize=>[$singleExonStep/10,1],
		   min=>[0,0],
		   max=>[1,0],
		   type=>"distribution",
		   tempName=>"tmp-single.distr",
		   stableName=>"grape-single.distr",
		   configMember=>"single-exon-lengths"
		  });
      }
    return $state;
  }
#===============================================================
sub installDistr
  {
    my ($point,$state,$config,$attributeName)=@_;
    my $controlParms=getMember($state,$attributeName);
    my $distr=getMember($point,$attributeName);
    my $filename=$controlParms->{tempName};
    my $configMember=$controlParms->{configMember};
    my $path="$dir/$filename";
    $config->{$configMember}=$path;
    saveDistr($distr,$filename);
  }
#===============================================================
sub installPoint
  {
    my ($point,$config,$state)=@_;
    $config->{"mean-intron-length"}=
      getMember($point,"intronLength");
    $config->{"mean-intergenic-length"}=
      getMember($point,"intergenicLength");
    $config->{"mean-5'-UTR-length"}=
      $config->{"mean-3'-UTR-length"}=getMember($point,"utrLength");
    $config->{"exon-optimism"}=
      getMember($point,"optimism");

    my $transFilename="$dir/grape.trans";
    writeTransFile($point,$transFilename);
    $config->{"transition-probabilities"}=$transFilename;

    if($optimizeDurations)
      {
	installDistr($point,$state,$config,"initialExonLengths");
	installDistr($point,$state,$config,"internalExonLengths");
	installDistr($point,$state,$config,"finalExonLengths");
	installDistr($point,$state,$config,"singleExonLengths");
      }
  }
#===============================================================
sub hashPoint
  {
    my ($point)=@_;
    my $nDim=@$point;
    my $hashString="";
    for(my $i=0 ; $i<$nDim ; ++$i)
      {
	my $value=$point->[$i];
	$hashString.="|$value";
      }
    return $hashString;
  }
#===============================================================
sub evaluatePoint
  {
    my ($point,$state)=@_;
    my $hashString=hashPoint($point);
    my $cachedValue=$pointScoreCache->{$hashString};
    if(defined($cachedValue)) {return $cachedValue}
    my $newConfig={}; %$newConfig=%$config;
    installPoint($point,$newConfig,$state);
    my $filestem="tmp";
    writeConfigFile($newConfig,$filestem);
    my $newScore=evaluate($filestem);
    $pointScoreCache->{$hashString}=$newScore;
    return $newScore;
  }
#===============================================================
sub gradientAscent
  {
    print LOG "Performing gradient ascent...\n";

    # Initialize some control structures
    $axisMap=initAxisMap();
    my $state=initState();

    # Initialize the current "point" in parameter space
    my $currentPoint=[];
    setMember($currentPoint,"intronLength",
	      $config->{"mean-intron-length"});
    setMember($currentPoint,"intergenicLength",
	      $config->{"mean-intergenic-length"});
    setMember($currentPoint,"utrLength",
	      $config->{"mean-5'-UTR-length"});
    setMember($currentPoint,"optimism",
	      $config->{"exon-optimism"});
    setMember($currentPoint,"ATG->GT",$transitions->{"ATG -> GT"});
    setMember($currentPoint,"AG->GT",$transitions->{"AG -> GT"});
    setMember($currentPoint,"polyA",$transitions->{"TAG -> POLYA"});
    setMember($currentPoint,"promoter",$transitions->{"-ATG -> -PROMOTER"});

    if($optimizeDurations)
      {
	setMember($currentPoint,"initialExonLengths",$initialExonLengths);
	setMember($currentPoint,"internalExonLengths",$internalExonLengths);
	setMember($currentPoint,"finalExonLengths",$finalExonLengths);
	setMember($currentPoint,"singleExonLengths",$singleExonLengths);
      }

    # Iterate until no changes have been made
    my $changes=1;
    my $currentScore=$baseline;
    $pointScoreCache->{hashPoint($currentPoint)}=$currentScore;
    my $numDim=@$currentPoint;
    while($changes)
      {
	$changes=0;

	# Iterate through all axes (dimensions)
	for(my $dim=0 ; $dim<$numDim ; ++$dim)
	  {
	    my $axisName=$axisMap->{$dim};
	    my $controlParms=$state->[$dim];
	    my $type=$controlParms->{type};
	    if($type eq "distribution")
	      {
		distributionStep(\$currentPoint,$state,$dim,
				 \$changes,\$currentScore);
		next;
	      }
	    if(!$controlParms->{direction})
	      {$controlParms->{direction}=(rand(1)<0.5?-1:1)}
	    my $min=$controlParms->{min}; my $max=$controlParms->{max};
	    my $delta=$controlParms->{direction}*$controlParms->{stepSize};
	    my $newPoint=[];@$newPoint=@$currentPoint;

	    # Try taking (another) step in the current direction...
	    my $newCoord=$currentPoint->[$dim]+$delta;
	    if($newCoord<=$max && $newCoord>=$min)
	      {
		$newPoint->[$dim]=$newCoord;
		my $newScore=evaluatePoint($newPoint,$state);
		if(compareScores($newScore,$currentScore)>0)
		  {
		    $currentPoint=$newPoint;
		    $currentScore=$newScore;
		    $changes=1;

		    printScore("1) $axisName=$newPoint->[$dim] " .
		      "delta=$delta ",$currentScore);
		    system("cp tmp.cfg best.cfg");
		    next;
		  }
	      }

	    # That direction is no good; try the other direction
	    $newCoord=$currentPoint->[$dim]-$delta;
	    if($newCoord<=$max && $newCoord>=$min)
	      {
		$newPoint->[$dim]=$newCoord;
		my $newScore=evaluatePoint($newPoint,$state);
		if(compareScores($newScore,$currentScore)>0)
		  {
		    $currentPoint=$newPoint;
		    $currentScore=$newScore;
		    $changes=1;
		    $controlParms->{direction}*=-1;

		    printScore("2) $axisName=$newPoint->[$dim] " .
		      "delta=$delta ",$currentScore);
		    system("cp tmp.cfg best.cfg");
		    next;
		  }
	      }

	    # Neither direction works, so reduce the step size & move on
	    # to the next axis (we'll come back to this axis next iteration)
	    $controlParms->{direction}=0;
	    if($controlParms->{stepSize}>$controlParms->{minStepSize})
	      {
		$controlParms->{stepSize}/=2;
		if($type eq "int")
		  {$controlParms->{stepSize}=int($controlParms->{stepSize})}
		if($controlParms->{stepSize}<$controlParms->{minStepSize})
		  {$controlParms->{stepSize}=$controlParms->{minStepSize}}
		$changes=1;

		print LOG "Reducing $axisName stepsize to " .
		  "$controlParms->{stepSize}\n";
	      }
	  }


	### check progress on the validation set

      }

    printScore("DONE.  Best score: ",$currentScore);
  }
#===============================================================
sub loadTransProbs
  {
    my ($filename)=@_;
    my $trans={};
    open(IN,$filename) || die "Can't open $filename\n";
    while(<IN>)
      {
	if(/^(\S+)\s*->\s*(\S+)\s*:\s*([\d\.]+)\s*$/)
	  {
	    my ($fromType,$toType,$P)=($1,$2,$3);
	    $trans->{"$fromType \-\> $toType"}=$P;
	  }
      }
    close(IN);
    return $trans;
  }
#===============================================================
sub writeTransFile
  {
    my ($point,$filename)=@_;

    my $ATG_GT=getMember($point,"ATG->GT");
    my $AG_GT=getMember($point,"AG->GT");
    my $polyA=getMember($point,"polyA");
    my $promoter=getMember($point,"promoter");

    my $ATG_TAG=1-$ATG_GT;
    my $AG_TAG=1-$AG_GT;
    my $nonPolyA=1-$polyA;
    my $nonPromoter=1-$promoter;
    my $TAG_PROM=0.5*$nonPolyA*$promoter;
    my $TAG_ATG=0.5*$nonPolyA*$nonPromoter;
    my $TAG_NEG_POLYA=0.5*$nonPolyA*$polyA;
    my $TAG_NEG_TAG=0.5*$nonPolyA*$nonPolyA;
    my $POLYA_PROM=0.5*$promoter;
    my $POLYA_ATG=0.5*$nonPromoter;
    my $POLYA_NEG_TAG=0.5*$nonPolyA;
    my $POLYA_NEG_POLYA=0.5*$polyA;
    my $NEG_ATG_ATG=0.5*$nonPromoter*$nonPromoter;
    my $NEG_ATG_PROM=0.5*$nonPromoter*$promoter;
    my $NEG_ATG_NEG_POLYA=0.5*$nonPromoter*$polyA;
    my $NEG_ATG_NEG_TAG=0.5*$nonPromoter*$nonPolyA;
    my $NEG_PROM_PROM=0.5*$promoter;
    my $NEG_PROM_ATG=0.5*$nonPromoter;
    my $NEG_PROM_NEG_TAG=0.5*$nonPolyA;
    my $NEG_PROM_NEG_POLYA=0.5*$polyA;

    my $ATG_GT0=$transitions->{"ATG -> GT0"};
    my $ATG_GT1=$transitions->{"ATG -> GT1"};
    my $ATG_GT2=$transitions->{"ATG -> GT2"};
    my $ATG_GT_SUM=$ATG_GT0+$ATG_GT1+$ATG_GT2;
    if($ATG_GT_SUM)
      {
	$ATG_GT0=$ATG_GT0/$ATG_GT_SUM*$ATG_GT;
	$ATG_GT1=$ATG_GT1/$ATG_GT_SUM*$ATG_GT;
	$ATG_GT2=$ATG_GT2/$ATG_GT_SUM*$ATG_GT;
      }
    else {$ATG_GT0=$ATG_GT1=$ATG_GT2=1/3}

    my $AG0_GT0=$transitions->{"AG0 -> GT0"};
    my $AG0_GT1=$transitions->{"AG0 -> GT1"};
    my $AG0_GT2=$transitions->{"AG0 -> GT2"};
    my $AG0_GT_SUM=$AG0_GT0+$AG0_GT1+$AG0_GT2;
    if($AG0_GT_SUM)
      {
	$AG0_GT0=$AG0_GT0/$AG0_GT_SUM*$AG_GT;
	$AG0_GT1=$AG0_GT1/$AG0_GT_SUM*$AG_GT;
	$AG0_GT2=$AG0_GT2/$AG0_GT_SUM*$AG_GT;
      }
    else {$AG0_GT0=$AG0_GT1=$AG0_GT2=1/3}

    my $AG1_GT0=$transitions->{"AG1 -> GT0"};
    my $AG1_GT1=$transitions->{"AG1 -> GT1"};
    my $AG1_GT2=$transitions->{"AG1 -> GT2"};
    my $AG1_GT_SUM=$AG1_GT0+$AG1_GT1+$AG1_GT2;
    if($AG1_GT_SUM)
      {
	$AG1_GT0=$AG1_GT0/$AG1_GT_SUM*$AG_GT;
	$AG1_GT1=$AG1_GT1/$AG1_GT_SUM*$AG_GT;
	$AG1_GT2=$AG1_GT2/$AG1_GT_SUM*$AG_GT;
      }
    else {$AG1_GT0=$AG1_GT1=$AG1_GT2=1/3}

    my $AG2_GT0=$transitions->{"AG2 -> GT0"};
    my $AG2_GT1=$transitions->{"AG2 -> GT1"};
    my $AG2_GT2=$transitions->{"AG2 -> GT2"};
    my $AG2_GT_SUM=$AG2_GT0+$AG2_GT1+$AG2_GT2;
    if($AG2_GT_SUM)
      {
	$AG2_GT0=$AG2_GT0/$AG2_GT_SUM*$AG_GT;
	$AG2_GT1=$AG2_GT1/$AG2_GT_SUM*$AG_GT;
	$AG2_GT2=$AG2_GT2/$AG2_GT_SUM*$AG_GT;
      }
    else {$AG2_GT0=$AG2_GT1=$AG2_GT2=1/3}

    my $NEG_GT_AG=$AG_GT;
    my $NEG_GT_ATG=$AG_TAG;

    my $NEG_GT0_AG0=$transitions->{"-GT0 -> -AG0"};
    my $NEG_GT0_AG1=$transitions->{"-GT0 -> -AG1"};
    my $NEG_GT0_AG2=$transitions->{"-GT0 -> -AG2"};
    my $NEG_GT0_AG_SUM=$NEG_GT0_AG0+$NEG_GT0_AG1+$NEG_GT0_AG2;
    if($NEG_GT0_AG_SUM)
      {
	$NEG_GT0_AG0=$NEG_GT0_AG0 / $NEG_GT0_AG_SUM * $NEG_GT_AG;
	$NEG_GT0_AG1=$NEG_GT0_AG1 / $NEG_GT0_AG_SUM * $NEG_GT_AG;
	$NEG_GT0_AG2=$NEG_GT0_AG2 / $NEG_GT0_AG_SUM * $NEG_GT_AG;
      }
    else {$NEG_GT0_AG0=$NEG_GT0_AG1=$NEG_GT0_AG2=1/3}

    my $NEG_GT1_AG0=$transitions->{"-GT1 -> -AG0"};
    my $NEG_GT1_AG1=$transitions->{"-GT1 -> -AG1"};
    my $NEG_GT1_AG2=$transitions->{"-GT1 -> -AG2"};
    my $NEG_GT1_AG_SUM=$NEG_GT1_AG0+$NEG_GT1_AG1+$NEG_GT1_AG2;
    if($NEG_GT1_AG_SUM)
      {
	$NEG_GT1_AG0=$NEG_GT1_AG0 / $NEG_GT1_AG_SUM * $NEG_GT_AG;
	$NEG_GT1_AG1=$NEG_GT1_AG1 / $NEG_GT1_AG_SUM * $NEG_GT_AG;
	$NEG_GT1_AG2=$NEG_GT1_AG2 / $NEG_GT1_AG_SUM * $NEG_GT_AG;
      }
    else {$NEG_GT1_AG0=$NEG_GT1_AG1=$NEG_GT1_AG2=1/3}
    
    my $NEG_GT2_AG0=$transitions->{"-GT2 -> -AG0"};
    my $NEG_GT2_AG1=$transitions->{"-GT2 -> -AG1"};
    my $NEG_GT2_AG2=$transitions->{"-GT2 -> -AG2"};
    my $NEG_GT2_AG_SUM=$NEG_GT2_AG0+$NEG_GT2_AG1+$NEG_GT2_AG2;
    if($NEG_GT2_AG_SUM)
      {
	$NEG_GT2_AG0=$NEG_GT2_AG0 / $NEG_GT2_AG_SUM * $NEG_GT_AG;
	$NEG_GT2_AG1=$NEG_GT2_AG1 / $NEG_GT2_AG_SUM * $NEG_GT_AG;
	$NEG_GT2_AG2=$NEG_GT2_AG2 / $NEG_GT2_AG_SUM * $NEG_GT_AG;
      }
    else {$NEG_GT2_AG0=$NEG_GT2_AG1=$NEG_GT2_AG2=1/3}

    my $NEG_TAG_AG=$ATG_GT;
    my $NEG_TAG_ATG=$ATG_TAG;

    my $NEG_TAG_AG0=$transitions->{"-TAG -> -AG0"};
    my $NEG_TAG_AG1=$transitions->{"-TAG -> -AG1"};
    my $NEG_TAG_AG2=$transitions->{"-TAG -> -AG2"};
    my $NEG_TAG_AG_SUM=$NEG_TAG_AG0+$NEG_TAG_AG1+$NEG_TAG_AG2;
    if($NEG_TAG_AG_SUM)
      {
	$NEG_TAG_AG0=$NEG_TAG_AG0 / $NEG_TAG_AG_SUM * $NEG_TAG_AG;
	$NEG_TAG_AG1=$NEG_TAG_AG1 / $NEG_TAG_AG_SUM * $NEG_TAG_AG;
	$NEG_TAG_AG2=$NEG_TAG_AG2 / $NEG_TAG_AG_SUM * $NEG_TAG_AG;
      }
    else {$NEG_TAG_AG0=$NEG_TAG_AG1=$NEG_TAG_AG2=1/3}

    open(OUT,">$filename") || die "can't write to file: $filename\n";
    print OUT "24\n"; # number of signal types
    print OUT "
PROMOTER -> ATG : 1
ATG -> GT : $ATG_GT
ATG -> TAG : $ATG_TAG
GT -> AG : 1
AG -> GT : $AG_GT
AG -> TAG : $AG_TAG
TAG -> POLYA : $polyA
TAG -> ATG : $TAG_ATG
TAG -> PROMOTER : $TAG_PROM
TAG -> -TAG : $TAG_NEG_TAG
TAG -> -POLYA : $TAG_NEG_POLYA
POLYA -> ATG : $POLYA_ATG
POLYA -> PROMOTER : $POLYA_PROM
POLYA -> -TAG : $POLYA_NEG_TAG
POLYA -> -POLYA : $POLYA_NEG_POLYA

ATG -> GT0 : $ATG_GT0
ATG -> GT1 : $ATG_GT1
ATG -> GT2 : $ATG_GT2
GT0 -> AG0 : 1
GT1 -> AG1 : 1
GT2 -> AG2 : 1

AG0 -> GT0 : $AG0_GT0
AG0 -> GT1 : $AG0_GT1
AG0 -> GT2 : $AG0_GT2
AG0 -> TAG : $AG_TAG

AG1 -> GT0 : $AG1_GT0
AG1 -> GT1 : $AG1_GT1
AG1 -> GT2 : $AG1_GT2
AG1 -> TAG : $AG_TAG

AG2 -> GT0 : $AG2_GT0
AG2 -> GT1 : $AG2_GT1
AG2 -> GT2 : $AG2_GT2
AG2 -> TAG : $AG_TAG

-PROMOTER -> ATG : $NEG_PROM_ATG
-PROMOTER -> PROMOTER : $NEG_PROM_PROM
-PROMOTER -> -TAG : $NEG_PROM_NEG_TAG
-PROMOTER -> -POLYA : $NEG_PROM_NEG_POLYA
-ATG -> -PROMOTER : $promoter
-ATG -> ATG : $NEG_ATG_ATG
-ATG -> PROMOTER : $NEG_ATG_PROM
-ATG -> -TAG : $NEG_ATG_NEG_TAG
-ATG -> -POLYA : $NEG_ATG_NEG_POLYA
-GT -> -ATG : $AG_TAG
-GT -> -AG : $AG_GT
-AG -> -GT : 1
-TAG -> -ATG : $NEG_TAG_ATG
-TAG -> -AG : $NEG_TAG_AG
-POLYA -> -TAG : 1

-AG0 -> -GT0 : 1
-AG1 -> -GT1 : 1
-AG2 -> -GT2 : 1

-GT0 -> -AG0 : $NEG_GT0_AG0
-GT0 -> -AG1 : $NEG_GT0_AG1
-GT0 -> -AG2 : $NEG_GT0_AG2
-GT0 -> -ATG : $NEG_GT_ATG

-GT1 -> -AG0 : $NEG_GT1_AG0
-GT1 -> -AG1 : $NEG_GT1_AG1
-GT1 -> -AG2 : $NEG_GT1_AG2
-GT1 -> -ATG : $NEG_GT_ATG

-GT2 -> -AG0 : $NEG_GT2_AG0
-GT2 -> -AG1 : $NEG_GT2_AG1
-GT2 -> -AG2 : $NEG_GT2_AG2
-GT2 -> -ATG : $NEG_GT_ATG

-TAG -> -AG0 : $NEG_TAG_AG0
-TAG -> -AG1 : $NEG_TAG_AG1
-TAG -> -AG2 : $NEG_TAG_AG2
";

    close(OUT);
  }
#===============================================================
# A distribution is an array of 2-element arrays, [x,y]
sub loadDistr
  {
    my ($filename)=@_;
    my $distr=[];
    open(IN,$filename);
    while(<IN>)
      {
	if(/^(\d+)\s+(\S+)/)
	  {
	    my ($x,$y)=($1,$2);
	    push(@$distr,[$x,$y]);
	  }
      }
    close(IN);
    return $distr;
  }
#===============================================================
# NOTE: |factor| should be between 0 and 1 (ideally around 5% or so)
sub adjustKurtosis
  {
    my ($originalDistr,$factor)=@_;
    my $distr=cloneDistr($originalDistr);
    my $n=@$distr;
    my $sum;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $elem=$distr->[$i];
	my $orig=$elem->[1];
	my $new=$orig+$factor;
	if($new<0) {$new=0}
	if($new<1) {$elem->[1]=$new}
	$sum+=$elem->[1];
      }
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $elem=$distr->[$i];
	$elem->[1]/=$sum;
      }
    return $distr;
  }
#===============================================================
sub adjustSkew
  {
    my ($originalDistr,$delta)=@_;
    my $distr=cloneDistr($originalDistr);
    my $n=@$originalDistr;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $fromIndex=$i-$delta;
	if($fromIndex<0)
	  {$distr->[$i]->[1]=$originalDistr->[0]->[1]/($delta-$i+1)}
	elsif($fromIndex>=$n)
	  {$distr->[$i]->[1]=$originalDistr->[$n-1]->[1]/
	     2**($fromIndex-$n+1)}
	else
	  {$distr->[$i]->[1]=$originalDistr->[$fromIndex]->[1]}
      }
    my $sum;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $elem=$distr->[$i];
	$sum+=$elem->[1];
      }
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $elem=$distr->[$i];
	$elem->[1]/=$sum;
      }
    return $distr;
  }
#===============================================================
sub saveDistr
  {
    my ($distr,$filename)=@_;
    open(OUT,">$filename") || die "Can't write file: $filename\n";
    my $n=@$distr;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $elem=$distr->[$i];
	my ($x,$y)=@$elem;
	print OUT "$x\t$y\n";
      }
    close(OUT);
  }
#===============================================================
sub cloneDistr
  {
    my ($distr)=@_;
    my $clone=[];
    my $n=@$distr;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $elem=$distr->[$i];
	my $newElem=[];
	@$newElem=@$elem;
	push @$clone,$newElem;
      }
    return $clone;
  }
#===============================================================
sub distributionStep
  {
    my ($currentPoint,$state,$dim,$changes,$currentScore)=@_;

    if(!kurtosisStep($currentPoint,$state,$dim,$changes,$currentScore))
      {skewStep($currentPoint,$state,$dim,$changes,$currentScore)}
  }
#===============================================================
sub kurtosisStep
  {
    my ($currentPoint,$state,$dim,$changes,$currentScore)=@_;

    my $controlParms=$state->[$dim];
    my $axisName=$axisMap->{$dim};
    if(!$controlParms->{direction}->[0])
      {$controlParms->{direction}->[0]=(rand(1)<0.5?-1:1)}
    my $stepSize=$controlParms->{stepSize}->[0];
    my $delta=$controlParms->{direction}->[0]*$stepSize;
    my $newPoint=[];@$newPoint=@$$currentPoint;
    my $distr=$$currentPoint->[$dim];
    my $tempFilename=$controlParms->{tempName};
    my $stableFilename=$controlParms->{stableName};
    my $configMember=$controlParms->{configMember};

    # Try taking a step in the current direction...
    my $newDistr=adjustKurtosis($distr,$delta);
    $newPoint->[$dim]=$newDistr;
    my $newScore=evaluatePoint($newPoint,$state);
    if(compareScores($newScore,$$currentScore)>0)
      {
	$$currentPoint=$newPoint;
	$$currentScore=$newScore;
	$$changes=1;
	system("mv $tempFilename $stableFilename");
	$config->{$configMember}="$dir/$stableFilename";
	printScore("$axisName kurtosis += $delta ",$$currentScore);
	return 1;
      }

    # That direction is no good; try the other direction
    my $newDistr=adjustKurtosis($distr,-$delta);
    $newPoint->[$dim]=$newDistr;
    my $newScore=evaluatePoint($newPoint,$state);
    if(compareScores($newScore,$$currentScore)>0)
      {
	$$currentPoint=$newPoint;
	$$currentScore=$newScore;
	$$changes=1;
	system("mv $tempFilename $stableFilename");
	$config->{$configMember}="$dir/$stableFilename";
	$controlParms->{direction}->[0]*=-1;
	printScore("$axisName kurtosis -= $delta ",$$currentScore);
	return 1;
      }

    # Neither direction works, so reduce the step size & move on
    # to the next axis (we'll come back to this axis next iteration)
    my $minStepSize=$controlParms->{minStepSize}->[0];
    if($stepSize>$minStepSize)
      {
	$stepSize/=2;
	if($stepSize<$minStepSize) {$stepSize=$minStepSize}
	$$changes=1;
	$controlParms->{stepSize}->[0]=$stepSize;
	$controlParms->{direction}->[0]=0;
	print LOG "Reducing $axisName kurtosis stepsize to $stepSize\n";
      }
    return 0;
  }
#===============================================================
sub skewStep
  {
    my ($currentPoint,$state,$dim,$changes,$currentScore)=@_;

    my $controlParms=$state->[$dim];
    my $axisName=$axisMap->{$dim};
    if(!$controlParms->{direction}->[1])
      {$controlParms->{direction}->[1]=(rand(1)<0.5?-1:1)}
    my $stepSize=$controlParms->{stepSize}->[1];
    my $delta=$controlParms->{direction}->[1]*$stepSize;
    my $newPoint=[];@$newPoint=@$$currentPoint;
    my $distr=$$currentPoint->[$dim];
    my $tempFilename=$controlParms->{tempName};
    my $stableFilename=$controlParms->{stableName};
    my $configMember=$controlParms->{configMember};

    # Try taking a step in the current direction...
    my $newDistr=adjustSkew($distr,$delta);
    $newPoint->[$dim]=$newDistr;
    my $newScore=evaluatePoint($newPoint,$state);
    if(compareScores($newScore,$$currentScore)>0)
      {
	$$currentPoint=$newPoint;
	$$currentScore=$newScore;
	$$changes=1;
	system("mv $tempFilename $stableFilename");
	$config->{$configMember}="$dir/$stableFilename";
	printScore("$axisName skew += $delta ",$$currentScore);
	return;
      }

    # That direction is no good; try the other direction
    my $newDistr=adjustSkew($distr,-$delta);
    $newPoint->[$dim]=$newDistr;
    my $newScore=evaluatePoint($newPoint,$state);
    if(compareScores($newScore,$$currentScore)>0)
      {
	$$currentPoint=$newPoint;
	$$currentScore=$newScore;
	$$changes=1;
	system("mv $tempFilename $stableFilename");
	$config->{$configMember}="$dir/$stableFilename";
	$controlParms->{direction}->[1]*=-1;
	printScore("$axisName skew -= $delta ",$$currentScore);
	return;
      }

    # Neither direction works, so reduce the step size & move on
    # to the next axis (we'll come back to this axis next iteration)
    my $minStepSize=$controlParms->{minStepSize}->[1];
    if($stepSize>$minStepSize)
      {
	$stepSize/=2;
	$stepSize=int($stepSize+5/9);
	if($stepSize<$minStepSize) {$stepSize=$minStepSize}
	$$changes=1;
	$controlParms->{stepSize}->[1]=$stepSize;
	$controlParms->{direction}->[1]=0;
	print LOG "Reducing $axisName skew stepsize to $stepSize\n";
      }
  }
#===============================================================
sub min
  {
    my ($a,$b)=@_;
    return ($a<$b ? $a : $b);
  }
#===============================================================
sub determineIsochore
  {
    my $ls=`ls iso*-*.gff`;
    my @files=split/\s+/,$ls;
    my $numFiles=@files;
    if($numFiles>1) 
      {die "Must be at most one iso*-*.gff file in current directory!\n"}
    my $file=$files[0];
    $file=~/iso(\d+)-(\d+)\.gff/
      || die "Can't find iso*-*.gff file in current directory!\n";
    my $isochore="$1\-$2";
    return $isochore;
  }
#===============================================================
#===============================================================
#===============================================================






