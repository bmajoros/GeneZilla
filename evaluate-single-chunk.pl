#!/usr/bin/perl
use strict;
use lib('/home/bmajoros/genomics/perl');
use GffTranscriptReader;
use FastaReader;

my $genezilla="out";
my $cdna="chunks";
my $gffReader=new GffTranscriptReader;
my $genezillaFiles=getFiles($genezilla);
my $cdnaFiles=getFiles($cdna);
my $numGeneZilla=@$genezillaFiles;
my $numCdna=@$cdnaFiles;
my $totalBases;
my $totalCoding;
my @globalPredictions;
my @globalCdnas;
my $out;
my $seqLen;

my $numFiles=min($numGeneZilla,$numCdna);

open(PIPE,"grep source-version $genezilla/* |") || die;
my $line=<PIPE>;
close(PIPE);
$line=~/source-version\s*(\S+)/;
my $programName=$1;
pad(\$programName,9);

print STDOUT "#files=$numFiles\n\n";
print STDOUT "         |NUCLEOTIDE--|SPLICE-|ATG/TAG|EXONS------|EXACT\n";
print STDOUT "         |SN  SP   F  |SN  SP |SN  SP |SN  SP  F  |GENES\n";

evaluateExons($programName,$genezillaFiles,$numFiles,$cdnaFiles,$genezilla,
	      $cdna);
evaluateNuc($programName,$genezillaFiles,$numFiles,$cdnaFiles,$genezilla,
	    $cdna);
evaluateGenes($programName,$genezillaFiles,$numFiles,$cdnaFiles,$genezilla,
	      $cdna);

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
	my $fasta="chunks/$id.fasta";
	$seqLen+=FastaReader::getGenomicSize($fasta);
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
#----------------------------------------------------------------
sub evaluateNuc
  {
    my ($label)=@_;
    my ($totalTP, $totalTN, $totalFP, $totalFN);

    my $cdnaVector=makeVector(\@globalCdnas);
    my $predVector=makeVector(\@globalPredictions);

    my ($percentIdentity,$numErrors)=
      getPercentIdentity($predVector,$cdnaVector);
    my ($sensitivity,$specificity,$TP,$TN,$FP,$FN)=
      getSensAndSpecificity($predVector,$cdnaVector);

    $totalBases+=($TP+$TN+$FP+$FN);
    $totalCoding+=($TP+$FN);

    my $sensitivity=int(100*$TP/($TP+$FN)+5/9);
    my $specificity=int(100*$TP/($TP+$FP)+5/9);
    #my $accuracy=int(100*($TP+$TN)/($TP+$TN+$FP+$FN)+5/9);
    my $F=int(2*$sensitivity*$specificity/
	      ($sensitivity+$specificity)+5/9);

    pad(\$sensitivity,2);
    pad(\$specificity,3);
    pad(\$F,2);
    print STDOUT "$label|$sensitivity\% $specificity\% $F\%";
    print STDOUT "$out";
  }
#----------------------------------------------------------------
sub countIdenticalExons
  {
    my ($list1,$list2)=@_;
    my $identical;
    my $n1=@$list1;
    my $n2=@$list2;
    my $firstJ=0;
  OUTER:
    for(my $i=0 ; $i<$n1 ; ++$i)
      {
	my $exon1=$list1->[$i];
	next if $exon1->{matched};###
	my $begin1=$exon1->getBegin();
	my $end1=$exon1->getEnd();
      INNER:
	for(my $j=$firstJ ; $j<$n2 ; ++$j)
	  {
	    my $exon2=$list2->[$j];
	    my $begin2=$exon2->getBegin();
	    if($begin2>$end1) {next OUTER}
	    my $end2=$exon2->getEnd();
	    if($end2<$begin1) {next INNER}
	    next if $exon2->{matched};###
	    if($begin1==$begin2 && $end1==$end2)
	      {
		++$identical;
		$exon1->{matched}=1;
		$exon2->{matched}=1;
		$firstJ=$j+1;
	      }
	  }
      }
    my $outOf=@$list1;
    my $outOf2=@$list2;
    return ($identical,$outOf,$outOf2);
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
    @list=sort {$a->getBegin() <=> $b->getBegin()} @list;
    return \@list;
  }
#----------------------------------------------------------------
sub trimPredictions
  {
    # This function drops any prediction which does not overlap
    # a known gene -- we simply have no way of knowing whether the
    # prediction is right or wrong.

    my $n1=@globalPredictions;
    my $n2=@globalCdnas;
    my @keep;
    my $firstJ=0;
    for(my $i=0 ; $i<$n1 ; ++$i)
      {
	my $prediction=$globalPredictions[$i];
	my $predBegin=$prediction->getBegin();
	my $predEnd=$prediction->getEnd();
	for(my $j=$firstJ ; $j<$n2 ; ++$j)
	  {
	    my $cdna=$globalCdnas[$j];
	    my $cdnaBegin=$cdna->getBegin();
	    my $cdnaEnd=$cdna->getEnd();
	    if($cdnaBegin>$predEnd) {last}
	    elsif($cdnaEnd<$predBegin) {next}
	    push @keep,$prediction;
	    $firstJ=$j;
	    last;
	  }
      }
    @globalPredictions=@keep;
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
	next unless $predictions;
	push @globalPredictions,@$predictions;

	$file=~/(\d+)/ || next;
	my $cdnaFile="$1.gff";
	next unless -e "$cdnaDir/$cdnaFile" && -e "$programDir/$file";
	my $cDNAs=$gffReader->loadGFF("$cdnaDir/$cdnaFile");
	push @globalCdnas,@$cDNAs;
      }
    trimPredictions();

    my $cdnaExons=getExonList(\@globalCdnas);
    my $cdnaSplices=getSpliceSites($cdnaExons);
    my $cdnaCodons=getStartAndStopCodons($cdnaExons);

    my $predExons=getExonList(\@globalPredictions);
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

    my $spliceSens=($spliceTP+$spliceFN ? int(100*$spliceTP/($spliceTP+$spliceFN)+5/9) : 1);
    my $spliceSpec=($spliceTP+$spliceFP ? int(100*$spliceTP/($spliceTP+$spliceFP)+5/9) : 1);

    my $codonSens=int(100*$codonTP/($codonTP+$codonFN)+5/9);
    my $codonSpec=int(100*$codonTP/($codonTP+$codonFP)+5/9);

    my $exonSpecificity=int(100*$totalIdentical/$totalExons+5/9);
    my $exonSensitivity=int(100*$totalIdentical/$totalExonsCdna+5/9);
    my $F=$exonSensitivity+$exonSpecificity>0 ?
      int(2*$exonSensitivity*$exonSpecificity/
      ($exonSensitivity+$exonSpecificity)+5/9) : 0;
    pad(\$spliceSens,2);
    pad(\$spliceSpec,2);
    pad(\$codonSens,2);
    pad(\$codonSpec,2);
    pad(\$exonSensitivity,2);
    pad(\$exonSpecificity,2);
    pad(\$F,2);
    $out.="|$spliceSens\% $spliceSpec\%|$codonSens\% $codonSpec\%|$exonSensitivity\% $exonSpecificity\% $F\%";
  }
#----------------------------------------------------------------
sub evaluateGenes
  {
    my ($label,$files,$numFiles,$cdnaFiles,$programDir,$cdnaDir)=@_;
    my $totalIdentical=0;
    my $numCdnaGenes=@globalCdnas;
    NEXTGENE:
      for(my $ic=0 ; $ic<$numCdnaGenes ; ++$ic)
	{
	  my $cdnaGene=$globalCdnas[$ic];
	  my $cdnaExons=$cdnaGene->{exons};
	  my $numExons=@$cdnaExons;
	  for(my $i=0 ; $i<$numExons ; ++$i) 
	    {
	      if(!$cdnaExons->[$i]->{matched}) {next NEXTGENE}
	    }
	  ++$totalIdentical;
	}
    my $percent=int(100*$totalIdentical/$numCdnaGenes+5/9);
    pad(\$percent,2);
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
	unless($exonType=~
              /initial-exon|internal-exon|final-exon|single-exon/)
	  {die "UNKNOWN EXON TYPE: $exonType"}
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
sub pad
  {
    my ($ptr,$len)=@_;
    my $l=length $$ptr;
    my $extra=$len-$l;
    my $padding=' 'x$extra;
    $$ptr.=$padding;
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

#-----------------------------------------------------------------
sub max
  {
    my ($a,$b)=@_;
    return ($a>$b ? $a : $b);
  }
#-----------------------------------------------------------------
sub getPercentIdentity
  {
    my ($vector1,$vector2)=@_;
    my $identity;
    my $length;
    for(my $i=0 ; $i<$seqLen ; ++$i)
      {
	my $elem1=$vector1->[$i];
	my $elem2=$vector2->[$i];
	next unless defined($elem1) || defined($elem2);
	if($elem1==$elem2){++$identity}
	++$length;
      }
    return ($identity/$length,$length-$identity);
  }
#-----------------------------------------------------------------
sub getSensAndSpecificity
  {
    my ($hypoth,$correct)=@_;
    my ($TP,$TN,$FP,$FN)=(0,0,0,0);
    for(my $i=0 ; $i<$seqLen ; ++$i)
      {
	my $hypothElem=$hypoth->[$i];
	my $correctElem=$correct->[$i];
	next unless defined($hypothElem) || defined($correctElem);
	if($hypothElem==$correctElem)
	  {
	    if($correctElem==1){++$TP}
	    else {++$TN}
	  }
	else
	  {
	    if($correctElem==1){++$FN}
	    else {++$FP}
	  }
      }
    my $sensitivity=($TP+$FN>0 ? $TP/($TP+$FN) : 1);
    my $specificity=($TP+$FP>0 ? $TP/($TP+$FP) : 1);
    return ($sensitivity,$specificity,$TP,$TN,$FP,$FN);
  }
#-----------------------------------------------------------------
sub makeVector
  {
    my ($transcripts)=@_;
    my $vector=[];
    my $n=@$transcripts;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $transcript=$transcripts->[$i];
	my $exons=$transcript->{exons};
	my $n=@$exons;
	my $prevEnd;
	for(my $i=0 ; $i<$n ; ++$i)
	  {
	    my $exon=$exons->[$i];
	    my $begin=$exon->getBegin();
	    my $end=$exon->getEnd();
	    if(defined($prevEnd))
	      {
		for(my $j=$prevEnd ; $j<$begin ; ++$j)
		  {
		    $vector->[$j]=0;
		  }
	      }
	    for(my $j=$begin ; $j<$end ; ++$j)
	      {
		$vector->[$j]=1;
	      }
	    $prevEnd=$end;
	  }
      }
    return $vector;
  }
#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------
