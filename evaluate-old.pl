#!/usr/bin/perl
use strict;
use lib('/home/bmajoros/genomics/perl');
use GffTranscriptReader;

my $genezilla="out";
my $cdna="chunks";
my $gffReader=new GffTranscriptReader;
my $genezillaFiles=getFiles($genezilla);
my $cdnaFiles=getFiles($cdna);
my $numGeneZilla=@$genezillaFiles;
my $numCdna=@$cdnaFiles;
my $totalBases;
my $totalCoding;

my $numFiles=min($numGeneZilla,$numCdna);

print STDERR "#files=$numFiles\n\n";
print STDERR "\t\tnuc\tnuc\tnuc\tsplice\tsplice\tATG/TAG\tATG/TAG\texon\texon\texact\n";
print STDERR "\t\tsens\tspec\taccur\tsens\tspec\tsens\tspec\tsens\tspec\tgenes\n";

evaluateNuc("GeneZilla ",$genezillaFiles,$numFiles,$cdnaFiles,$genezilla,$cdna);
evaluateExons("GeneZilla ",$genezillaFiles,$numFiles,$cdnaFiles,$genezilla,$cdna);
evaluateGenes("GeneZilla ",$genezillaFiles,$numFiles,$cdnaFiles,$genezilla,$cdna);

#my $d=0.61;
#my $density=$totalCoding/$totalBases;
#my $baseline=int(100*$d*$density+(1-$d)*(1-$density)+0.5);
#print STDERR "\nbaseline nuc accuracy: $baseline\%\n";

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
	#my $cdnaFile=$cdnaFiles->[$i];
	$file=~/(\d+)/ || next;
	my $cdnaFile="$1.gff";
	my $substrateLength=getSubstrateLength($cdnaDir,$cdnaFile);
	my $command="cmp-gff.pl $programDir/$file $cdnaDir/$cdnaFile $substrateLength";
	my $results=`$command`;
	$results=~/TP=(\d+)\s+TN=(\d+)\s+FP=(\d+)\s+FN=(\d+)/ || die $results;
	my ($TP,$TN,$FP,$FN)=($1,$2,$3,$4);

	$totalTP+=$TP;
	$totalTN+=$TN;
	$totalFP+=$FP;
	$totalFN+=$FN;
      }

    $totalBases+=($totalTP+$totalTN+$totalFP+$totalFN);
    $totalCoding+=($totalTP+$totalFN);

    my $sensitivity=int(100*$totalTP/($totalTP+$totalFN)+0.5);
    my $specificity=int(100*$totalTP/($totalTP+$totalFP)+0.5);
    my $accuracy=int(100*($totalTP+$totalTN)/
		     ($totalTP+$totalTN+$totalFP+$totalFN)+0.5);
    print STDERR "$label:\t$sensitivity\%\t$specificity\%\t$accuracy\%";
  }
#----------------------------------------------------------------
sub getSubstrateLength
  {
    my ($dir,$filename)=@_;
    my $transcripts=$gffReader->loadGFF("$dir/$filename");
    my $length=0;
    foreach my $transcript (@$transcripts)
      {
	my $m=max($transcript->getBegin(),$transcript->getEnd());
	if($m>$length) {$length=$m}
      }
    return $length;###+100;
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
	my $predExons=getExonList($predictions);
	my $predSplices=getSpliceSites($predExons);
	my $predCodons=getStartAndStopCodons($predExons);

	$file=~/(\d+)/ || next;
	my $cdnaFile="$1.gff";
	#my $cdnaFile=$cdnaFiles->[$i];
	my $cDNAs=$gffReader->loadGFF("$cdnaDir/$cdnaFile");
	my $cdnaExons=getExonList($cDNAs);
	my $cdnaSplices=getSpliceSites($cdnaExons);
	my $cdnaCodons=getStartAndStopCodons($cdnaExons);

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
	#print "splice: $commonSplices $spliceTP $spliceFP $spliceFN\n";

	my $commonCodons=
	  countCommonSpliceSites($predCodons,$cdnaCodons);
	$codonTP+=$commonCodons;
	$codonFP+=(@$predCodons-$commonCodons);
	$codonFN+=(@$cdnaCodons-$commonCodons);
      }

    my $spliceSens=int(100*$spliceTP/($spliceTP+$spliceFN)+0.5);
    my $spliceSpec=int(100*$spliceTP/($spliceTP+$spliceFP)+0.5);

    my $codonSens=int(100*$codonTP/($codonTP+$codonFN)+0.5);
    my $codonSpec=int(100*$codonTP/($codonTP+$codonFP)+0.5);

    my $exonSpecificity=int(100*$totalIdentical/$totalExons+0.5);
    my $exonSensitivity=int(100*$totalIdentical/$totalExonsCdna+0.5);
    print STDERR "\t$spliceSens\%\t$spliceSpec\%\t$codonSens\%\t$codonSpec\%\t$exonSensitivity\%\t$exonSpecificity\%";
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
	my $predGene=$predictions->[0];
	my $predExons=$predGene->{exons};

	$file=~/(\d+)/ || next;
	my $cdnaFile="$1.gff";
	#my $cdnaFile=$cdnaFiles->[$i];
	my $cDNAs=$gffReader->loadGFF("$cdnaDir/$cdnaFile");
	my $cdnaGene=$cDNAs->[0];
	my $cdnaExons=$cdnaGene->{exons};
	
	my ($identical,$outOf)=
	  countIdenticalExons($predExons,$cdnaExons);
	if($identical==$outOf && @$predictions==@$cDNAs &&
	   @$predExons==@$cdnaExons) 
	  {
	    ++$totalIdentical;
	    #print "\n$i EXACTGENE: $label ".(0+@$predExons)."\n";
	  }
	#else {print "\n$i WRONGGENE: $label $identical identical out of $outOf\n"}
      }
    my $percent=int(100*$totalIdentical/$totalGenes+0.5);
    print STDERR "\t$percent\%\n";
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

