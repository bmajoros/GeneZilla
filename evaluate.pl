#!/usr/bin/perl
use strict;
use GffTranscriptReader;
$|=1;

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

my $numFiles=$numCdna;#min($numGeneZilla,$numCdna);

print "Warning: This program works only if the test genes are in separate
chunk files -- if instead all the genes are in a single file, please use
evaluate-single-chunk.pl instead of this program.\n\n";

print STDOUT "#files=$numFiles\n";
print STDOUT "         |NUCLEOTIDE--|SPLICE-|ATG/TAG|EXONS------|EXACT\n";
print STDOUT "         |SN  SP   F  |SN  SP |SN  SP |SN  SP  F  |GENES\n";

evaluateNuc($programName,$genezillaFiles,$numFiles,$cdnaFiles,$genezilla,
	    $cdna);
evaluateExons($programName,$genezillaFiles,$numFiles,$cdnaFiles,$genezilla,
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

    my $spliceSens=$spliceTP+$spliceFN>0 ?
      int(100*$spliceTP/($spliceTP+$spliceFN)+5/9) : 0;
    my $spliceSpec=$spliceTP+$spliceFP>0 ?
      int(100*$spliceTP/($spliceTP+$spliceFP)+5/9) : 0;

    my $codonSens=$codonTP+$codonFN>0 ?
      int(100*$codonTP/($codonTP+$codonFN)+5/9) : 0;
    my $codonSpec=$codonTP+$codonFP>0 ?
      int(100*$codonTP/($codonTP+$codonFP)+5/9) : 0;

    my $exonSpecificity=int(100*$totalIdentical/$totalExons+5/9);
    my $exonSensitivity=int(100*$totalIdentical/$totalExonsCdna+5/9);
    my $F=$exonSensitivity+$exonSpecificity>0 ?
      int(2*$exonSensitivity*$exonSpecificity/
	  ($exonSensitivity+$exonSpecificity)+5/9) : 0;
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

