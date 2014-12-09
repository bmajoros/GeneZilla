#!/usr/bin/perl
######################################################################
# 4/29/2005 bmajoros@tigr.org
#
# New evaluation script which accounts for alternative transcripts
# and overlapping genes in the reference set (but not in the predict-
# ions).
######################################################################
use strict;
use GffTranscriptReader;

my $programName="GeneZilla";
if(@ARGV==1) {$programName=shift @ARGV}

my $genezilla="out";
my $cdna="chunks";
my $gffReader=new GffTranscriptReader;
my $genezillaFiles=getFiles($genezilla);
my $cdnaFiles=getFiles($cdna);
my $numGeneZilla=@$genezillaFiles;
my $numCdna=@$cdnaFiles;
my $totalBases;
my $totalCoding;
my ($nucSn,$nucSp,$nucF,$exonSn,$exonSp,$exonF);
my ($nTP,$nFP,$nTN,$nFN);
my ($eTP,$eFP,$eTN,$eFN);
my ($gTP,$gFP,$gTN,$gFN);

my $numFiles=$numCdna;

print STDOUT "#files=$numFiles\n";
print STDOUT "         |NUCLEOTIDE--|SPLICE-|ATG/TAG|EXONS------|EXACT\n";
print STDOUT "         |SN  SP   F  |SN  SP |SN  SP |SN  SP  F  |GENES\n";

#=====================================================================
# Iterate through all files, collecting counts that will later be used
# to compute accuracy scores
#=====================================================================
for(my $i=0 ; $i<$numFiles ; ++$i)
  {
    my $file=$cdnaFiles->[$i];
    $file=~/(\d+)/ || next;
    my $referenceFile="chunks/$1.gff";
    my $predictionFile="out/$1.gff";
    print STDERR "processing $referenceFile vs. $predictionFile\n";

    #======================================
    # Get the lists of genes for this chunk
    #======================================
    my $refGenes=$gffReader->loadGenes($referenceFile);
    my $predTranscripts=$gffReader->loadGFF($predictionFile);

    #============================================================
    # For each predicted transcript, identify all those reference
    # genes which overlap it
    #============================================================
    my $refIndex=0;
    my $numRefGenes=@$refGenes;
    my $numPredTranscripts=@$predTranscripts;
    print "$numRefGenes refGenes, $numPredTranscripts predTranscripts\n";
    for(my $i=0 ; $i<$numPredTranscripts ; ++$i)
      {
	my $predTrans=$predTranscripts->[$i];
	my $predBegin=$predTrans->getBegin();
	my $predEnd=$predTrans->getEnd();
	print "==========> $predBegin,$predEnd\n";
	print "REFINDEX=$refIndex\n";
	my @overlappingRefGenes;
	for(my $j=$refIndex ; $j<$numRefGenes ; ++$j)
	  {
	    my $refGene=$refGenes->[$j];
	    print "XXX $predBegin,$predEnd vs. " . $refGene->getBegin() ."," . $refGene->getEnd() . "\n";
	    if($predBegin>=$refGene->getEnd()) {++$refIndex; next} # speed!
	    if($predEnd<$refGene->getBegin()) {last}
	    push @overlappingRefGenes,$refGene;
	  }

	#===================================================================
	# If any of the overlapping genes overlap each other, quit (for now)
	#===================================================================
	my $numOverlapping=@overlappingRefGenes;
	print "$numOverlapping OVERLAPPING\n";
	next unless $numOverlapping;
	for(my $i=0 ; $i<$numOverlapping ; ++$i)
	  {
	    my $geneI=$overlappingRefGenes[$i];
	    for(my $j=$i+1 ; $j<$numOverlapping ; ++$j)
	      {
		my $geneJ=$overlappingRefGenes[$j];
		if($geneI->getBegin()<$geneJ->getEnd() &&
		   $geneI->getEnd()>$geneJ->getBegin())
		  {die "overlapping reference genes not handled yet\n"}
	      }
	  }
	
	#================================================================
	# For each reference gene overlapping this predicted gene, choose
	# the "best" transcript (i.e., the one that matches best with the
	# predicted transcript)
	#================================================================
	my @referenceTranscripts;
	for(my $i=0 ; $i<$numOverlapping ; ++$i)
	  {
	    my $refGene=$overlappingRefGenes[$i];
	    my $numRefTranscripts=$refGene->getNumTranscripts();
	    if($numRefTranscripts==1) 
	      {push @referenceTranscripts,$refGene->getIthTranscript(0)}
	    else
	      {
		my ($bestTrans,$bestF);
		for(my $j=0 ; $j<$numRefTranscripts ; ++$j)
		  {
		    my $refTranscript=$refGene->getIthTranscript($i);
		    my ($F,$Sn,$Sp,$TP,$TN,$FP,$FN)=
		      getNucAccuracy($refTranscript,$predTrans);
		    if(!defined($bestF) || $F>$bestF)
		      {$bestTrans=$refTranscript; $bestF=$F}
		  }
		push @referenceTranscripts,$bestTrans;
	      }
	  }
	my $numRefs=@referenceTranscripts;
	print "$numRefs reference genes overlap this prediction\n";
	next unless $numRefs>0;

	#===========================================================
	# Compute TP/FP/TN/FN statistics at all levels between these
	# transcripts, and use these to update the global statistics
	#===========================================================

	# First, nucleotides:
	my ($absBegin,$absEnd)=
	  getAbsBeginAndEnd($predTrans,\@referenceTranscripts);
	my $refVector=
	  makeMultiRefVector($absBegin,$absEnd-$absBegin,\@referenceTranscripts);
	my $predVector=makeVector($absBegin,$absEnd-$absBegin,$predTrans);

	my $debug=@$refVector;
	print "$absBegin,$absEnd,$debug\n";

	my ($F,$Sn,$Sp,$TP,$TN,$FP,$FN)=
	  nucAccuracyFromVectors($refVector,$predVector,$absEnd);
	$nTP+=$TP;
	$nTN+=$TN;
	$nFP+=$FP;
	$nFN+=$FN;
	print "$TP $TN $FP $FN\n";

	# Now exons:


	# Finally, whole genes:



      }


  }

#==========================
# Compute global statistics
#==========================
my $sensitivity=int(100*$nTP/($nTP+$nFN)+5/9);
my $specificity=int(100*$nTP/($nTP+$nFP)+5/9);
my $F=int(2*$sensitivity*$specificity/($sensitivity+$specificity)+5/9);
print STDOUT "$programName|$sensitivity\% $specificity\%  $F\%";


#===============================
#evaluateNuc($programName,$genezillaFiles,$numFiles,$cdnaFiles,$genezilla,
#	    $cdna);
#evaluateExons($programName,$genezillaFiles,$numFiles,$cdnaFiles,$genezilla,
#	      $cdna);
#evaluateGenes($programName,$genezillaFiles,$numFiles,$cdnaFiles,$genezilla,
#	      $cdna);

#----------------------------------------------------------------
# my ($absBegin,$absEnd)=getAbsBeginAndEnd($predTrans,\@referenceTranscripts);
sub getAbsBeginAndEnd
  {
    my ($predTrans,$refTranses)=@_;
    my $n=@$refTranses;
    my $begin=$predTrans->getBegin();
    my $end=$predTrans->getEnd();
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $refTrans=$refTranses->[$i];
	my $rBegin=$refTrans->getBegin();
	my $rEnd=$refTrans->getEnd();
	if($rBegin<$begin) {$begin=$rBegin}
	if($rEnd>$end) {$end=$rEnd}
      }
    return ($begin,$end);
  }
#----------------------------------------------------------------
# my $refVector=makeMultiRefVector($absBegin,$len,\@referenceTranscripts);
sub makeMultiRefVector
  {
    my ($offset,$len,$transcripts)=@_;
    
    my $vector=[];
    for(my $i=0 ; $i<$len ; ++$i) {$vector->[$i]=0}
    my $numTranscripts=@$transcripts;
    for(my $i=0 ; $i<$numTranscripts ; ++$i)
      {
	my $transcript=$transcripts->[$i];
	my $numExons=$transcript->numExons();
	for(my $i=0 ; $i<$numExons ; ++$i)
	  {
	    my $exon=$transcript->getIthExon($i);
	    my $exonBegin=$exon->getBegin();
	    my $exonEnd=$exon->getEnd();
	    for(my $j=$exonBegin ; $j<$exonEnd ; ++$j)
	      {
		$vector->[$j-$offset]=1;
	      }
	  }
      }
    return $vector;
  }
#----------------------------------------------------------------
sub makeVector
  {
    my ($offset,$len,$transcript)=@_;
    my $vector=[];
    for(my $i=0 ; $i<$len ; ++$i) {$vector->[$i]=0}
    my $numExons=$transcript->numExons();
    for(my $i=0 ; $i<$numExons ; ++$i)
      {
	my $exon=$transcript->getIthExon($i);
	my $exonBegin=$exon->getBegin();
	my $exonEnd=$exon->getEnd();
	for(my $j=$exonBegin ; $j<$exonEnd ; ++$j)
	  {
	    $vector->[$j-$offset]=1;
	  }
      }
    return $vector;
  }
#----------------------------------------------------------------
#  my ($F,$Sn,$Sp,$TP,$TN,$FP,$FN)=
#                  nucAccuracyFromVectors($rVect,$pVect,$absEnd);
sub nucAccuracyFromVectors
  {
    my ($refVector,$predVector,$absEnd)=@_;
    my ($TP,$TN,$FP,$FN)=(0,0,0,0);
    for(my $i=0 ; $i<$absEnd ; ++$i)
      {
	my $refDigit=$refVector->[$i];
	my $predDigit=$predVector->[$i];
	if($refDigit==1)
	  {
	    if($predDigit==1) {++$TP}
	    else {++$FN}
	  }
	else
	  {
	    if($predDigit==1) {++$FP}
	    else {++$TN}
	  }
      }
    my $Sn=$TP/($TP+$FN);
    my $Sp=($TP+$FP>0 ? $TP/($TP+$FP) : 0);
    my $F=($Sn+$Sp>0 ? 2*$Sn*$Sp/($Sn+$Sp) : 0);
    return ($F,$Sn,$Sp,$TP,$TN,$FP,$FN);
  }
#----------------------------------------------------------------
#  my ($F,$Sn,$Sp,$TP,$TN,$FP,$FN)=getNucAccuracy($refTrans,$predTrans);
sub getNucAccuracy
  {
    my ($refTranscript,$predTranscript)=@_;
    my $rBegin=$refTranscript->getBegin();
    my $rEnd=$refTranscript->getEnd();
    my $pBegin=$predTranscript->getBegin();
    my $pEnd=$predTranscript->getEnd();
    my $absBegin=($rBegin<$pBegin ? $rBegin : $pBegin);
    my $absEnd=($rEnd>$pEnd ? $rEnd : $pEnd);
    my $absLen=$absEnd-$absBegin;
    my $refVector=makeVector($absBegin,$absLen,$refTranscript);
    my $predVector=makeVector($absBegin,$absLen,$predTranscript);
    return nucAccuracyFromVectors($refVector,$predVector,$absEnd);
  }
#----------------------------------------------------------------
sub getFiles
  {
    my ($dir)=@_;
    my $files=`ls $dir`;
    my @files=split/\s+/,$files;
    @files=sort {$a=~/(\d+)/;my $a1=$1;$b=~/(\d+)/;$a1<=>$b} @files;
    my $n=@files;
    my @list;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $file=$files[$i];
	next unless $file=~/(\d+)\.gff/;
	my $id=$1;
	push @list,$file;
      }
    return \@list;
  }
#----------------------------------------------------------------
sub min
  {
    my ($a,$b)=@_;
    return ($a<$b ? $a : $b);
  }
#----------------------------------------------------------------
sub max
  {
    my ($a,$b)=@_;
    return ($a>$b ? $a : $b);
  }










#################################################################
#                         OLD ROUTINES
#################################################################




#----------------------------------------------------------------
sub evaluateNuc
  {
    my ($label,$files,$numFiles,$cdnaFiles,$programDir,$cdnaDir)=@_;
    my ($totalTP, $totalTN, $totalFP, $totalFN);
    for(my $i=0 ; $i<$numFiles ; ++$i)
      {
	my $file=$files->[$i];
	$file=~/(\d+)/ || next;
	my $cdnaFile="$1.gff";
	my $command="cmp-gff.pl $programDir/$file $cdnaDir/$cdnaFile";
	my $results=`$command`;
	$results=~/TP=(\d+)\s+TN=(\d+)\s+FP=(\d+)\s+FN=(\d+)/ || next;
	my ($TP,$TN,$FP,$FN)=($1,$2,$3,$4);

	$totalTP+=$TP;
	$totalTN+=$TN;
	$totalFP+=$FP;
	$totalFN+=$FN;
      }

    $totalBases+=($totalTP+$totalTN+$totalFP+$totalFN);
    $totalCoding+=($totalTP+$totalFN);

    my $sensitivity=int(100*$totalTP/($totalTP+$totalFN)+5/9);
    my $specificity=int(100*$totalTP/($totalTP+$totalFP)+5/9);
    my $accuracy=int(100*($totalTP+$totalTN)/
		     ($totalTP+$totalTN+$totalFP+$totalFN)+5/9);
    my $F=int(2*$sensitivity*$specificity/($sensitivity+$specificity)+5/9);
    print STDOUT "$label|$sensitivity\% $specificity\%  $F\%";
  }
#----------------------------------------------------------------
sub countIdenticalExons
  {
    my ($list1,$list2)=@_;
    my $identical;
    foreach my $exon1 (@$list1)
      {
	foreach my $exon2 (@$list2)
	  {
	    if(exonsAreIdentical($exon1,$exon2)) {++$identical}
	  }
      }
    my $outOf=@$list1;
    my $outOf2=@$list2;
    return ($identical,$outOf,$outOf2);
  }
#----------------------------------------------------------------
sub exonsAreIdentical
  {
    my ($exon1,$exon2)=@_;
    return 
      $exon1->getBegin()==$exon2->getBegin() &&
      $exon1->getEnd()==$exon2->getEnd();
  }
#----------------------------------------------------------------
sub getExonList
  {
    my ($transcripts)=@_;
    my @list;
    foreach my $transcript (@$transcripts)
      {
	my $n=$transcript->numExons();
	for(my $i=0 ; $i<$n ; ++$i)
	  {
	    my $exon=$transcript->getIthExon($i);
	    push @list,$exon;
	  }
      }
    my $debug=@list;
    return \@list;
  }
#----------------------------------------------------------------
sub evaluateExons
  {
    my ($label,$files,$numFiles,$cdnaFiles,$programDir,$cdnaDir)=@_;
    my ($totalExons,$totalExonsCdna,$totalIdentical);
    my ($spliceTP,$spliceFP,$spliceFN);
    my ($codonTP,$codonFP,$codonFN);
    for(my $i=0 ; $i<$numFiles ; ++$i)
      {
	my $file=$files->[$i];
	my $predictions=$gffReader->loadGFF("$programDir/$file");

	$file=~/(\d+)/ || next;
	my $cdnaFile="$1.gff";
	my $cDNAs=$gffReader->loadGFF("$cdnaDir/$cdnaFile");

	my $cdnaExons=getExonList($cDNAs);
	my $cdnaSplices=getSpliceSites($cdnaExons);
	my $cdnaCodons=getStartAndStopCodons($cdnaExons);

	my $leftTerminus=getLeftTerminus($cDNAs);
	my $rightTerminus=getRightTerminus($cDNAs);

	dropSuperfluousFeatures($predictions,$leftTerminus,$rightTerminus);

	my $predExons=getExonList($predictions);
	my $predSplices=getSpliceSites($predExons);
	my $predCodons=getStartAndStopCodons($predExons);

	my ($identical,$outOf,$outOfCdna)=
	  countIdenticalExons($predExons,$cdnaExons);
	$totalIdentical+=$identical;
	$totalExons+=$outOf;
	$totalExonsCdna+=$outOfCdna;

	my $commonSplices=
	  countCommonSpliceSites($predSplices,$cdnaSplices);
	$spliceTP+=$commonSplices;
	$spliceFP+=(@$predSplices-$commonSplices);
	$spliceFN+=(@$cdnaSplices-$commonSplices);

	my $commonCodons=
	  countCommonSpliceSites($predCodons,$cdnaCodons);
	$codonTP+=$commonCodons;
	$codonFP+=(@$predCodons-$commonCodons);
	$codonFN+=(@$cdnaCodons-$commonCodons);
      }

    my $spliceSens=int(100*$spliceTP/($spliceTP+$spliceFN)+5/9);
    my $spliceSpec=int(100*$spliceTP/($spliceTP+$spliceFP)+5/9);

    my $codonSens=int(100*$codonTP/($codonTP+$codonFN)+5/9);
    my $codonSpec=int(100*$codonTP/($codonTP+$codonFP)+5/9);

    my $exonSpecificity=int(100*$totalIdentical/$totalExons+5/9);
    my $exonSensitivity=int(100*$totalIdentical/$totalExonsCdna+5/9);
    my $F=int(2*$exonSensitivity*$exonSpecificity/
	      ($exonSensitivity+$exonSpecificity)+5/9);
    print STDOUT "|$spliceSens\% $spliceSpec\%|$codonSens\% $codonSpec\%|$exonSensitivity\% $exonSpecificity\% $F\%";
  }
#----------------------------------------------------------------
sub evaluateGenes
  {
    my ($label,$files,$numFiles,$cdnaFiles,$programDir,$cdnaDir)=@_;
    my ($totalGenes,$totalIdentical);
    for(my $i=0 ; $i<$numFiles ; ++$i)
      {
	++$totalGenes;
	my $file=$files->[$i];
	my $predictions=$gffReader->loadGFF("$programDir/$file");
	next unless @$predictions>0;

	$file=~/(\d+)/ || next;
	my $cdnaFile="$1.gff";
	my $cDNAs=$gffReader->loadGFF("$cdnaDir/$cdnaFile");

	my $leftTerminus=getLeftTerminus($cDNAs);
	my $rightTerminus=getRightTerminus($cDNAs);

	my $cdnaGene=$cDNAs->[0];
	my $cdnaExons=$cdnaGene->{exons};
	
	dropSuperfluousFeatures($predictions,$leftTerminus,$rightTerminus);
	my $predGene=$predictions->[0];
	my $predExons=$predGene->{exons};

	my ($identical,$outOf)=
	  countIdenticalExons($predExons,$cdnaExons);
	if($identical==$outOf && @$predictions==@$cDNAs &&
	   @$predExons==@$cdnaExons) 
	  {
	    ++$totalIdentical;
	    #my $numExons=@$predExons;
	    #print "identical gene has $numExons exons\n";
	  }
      }
    $totalIdentical+=0; # so Perl knows it's a number, not a string
    my $percent=int(100*$totalIdentical/$totalGenes+5/9);
    print STDOUT "|$percent\% ($totalIdentical)\n";
  }
#----------------------------------------------------------------
sub getSpliceSites
  {
    my ($exons)=@_;
    my (@donors,@acceptors);
    foreach my $exon (@$exons)
      {
	my $exonType=$exon->getType();
	if($exonType eq "initial-exon" || $exonType eq "internal-exon")
	  {
	    push @donors,$exon->getEnd();
	  }
	if($exonType eq "final-exon" || $exonType eq "internal-exon")
	  {
	    push @acceptors,$exon->getBegin();
	  }
      }
    my @spliceSites;
    push @spliceSites,@donors;
    push @spliceSites,@acceptors;
    return \@spliceSites;
  }
#----------------------------------------------------------------
sub getStartAndStopCodons
  {
    my ($exons)=@_;
    my (@starts,@stops);
    foreach my $exon (@$exons)
      {
	my $exonType=$exon->getType();
	if($exonType eq "initial-exon" || $exonType eq "single-exon")
	  {
	    push @starts,$exon->getBegin();
	  }
	if($exonType eq "final-exon" || $exonType eq "single-exon")
	  {
	    push @stops,$exon->getEnd();
	  }
      }
    my @list;
    push @list,@starts;
    push @list,@stops;
    return \@list;
  }
#----------------------------------------------------------------
sub countCommonSpliceSites
  {
    my ($sites1,$sites2)=@_;
    my (%hash,$common);
    foreach my $site (@$sites1) {$hash{$site}=1}
    foreach my $site (@$sites2) {++$common if $hash{$site}}

    return $common;
  }
#-----------------------------------------------------------------
#my $leftTerminus=getLeftTerminus($correctFeatures);
sub getLeftTerminus
  {
    my ($features)=@_;
    my $left;
    foreach my $feature (@$features)
      {
	if(!defined($left) || $feature->getBegin()<$left)
	  { $left=$feature->getBegin() }
      }
    return $left;
  }
#-----------------------------------------------------------------
#my $rightTerminus=getRightTerminus($correctFeatures);
sub getRightTerminus
  {
    my ($features)=@_;
    my $right;
    foreach my $feature (@$features)
      {
	if(!defined($right) || $feature->getEnd()>$right)
	  { $right=$feature->getEnd() }
      }
    return $right;
  }
#----------------------------------------------------------------
#dropSuperfluousFeatures($predictions,$leftTerminus,$rightTerminus);
sub dropSuperfluousFeatures
  {
    my ($predictions,$leftBoundary,$rightBoundary)=@_;
    my @retval;
    my $n=@$predictions;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $feature=$predictions->[$i];
	my $begin=$feature->getBegin();
	my $end=$feature->getEnd();
	if($end<=$leftBoundary || $begin>=$rightBoundary) {next}
	push @retval,$feature;

	if(0)
	  {
	    my $ne=$feature->numExons();
	    for(my $i=0 ; $i<$ne ; ++$i)
	      {
		my $exon=$feature->getIthExon($i);
		my $begin=$exon->getBegin();
		my $end=$exon->getEnd();
		if($end<=$leftBoundary || $begin>=$rightBoundary)
		  {
		    #print "LT=$leftBoundary RT=$rightBoundary b=$begin e=$end\n";
		    $feature->deleteExon($i);
		    --$i;
		    --$ne;
		  }
	      }
	  }
      }
    @$predictions=@retval;
  }
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

