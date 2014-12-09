#!/usr/bin/perl
use strict;

##############################################################################
# NOTE: I used a uniform distribution for the non-consensus bases of the
# null model in determining a cutoff threshold -- this should probably
# be changed to infer a true background distribution from the other WMM 
# entries.
##############################################################################

my $NUM_STRINGS=500;

my $usage="$0 <transfac-matrix> <consensus-offset> <consensus-length> <min-sensitivity>\n\t<output-filestem>";
die "$usage\n" unless @ARGV==5;
my ($filename,$consensusOffset,$consensusLength,$minSensitivity,
    $filestem)=@ARGV;

my (@headers,%headers,@matrixOutput,@matrix,$nullWMM);
open(IN,$filename) || die;
while(<IN>)
  {
    if(/^\#/)
      {
	next;
      }
    elsif(/^[ATCG\s]+$/)
      {
	@headers=split/\s+/,$_;
	my $n=@headers;
	for(my $i=0 ; $i<$n ; ++$i) {$headers{$headers[$i]}=$i}
      }
    elsif(/\d/)
      {
	my @fields=split/\s+/,$_;
	my $A=safelog($fields[$headers{'A'}]);
	my $T=safelog($fields[$headers{'T'}]);
	my $C=safelog($fields[$headers{'C'}]);
	my $G=safelog($fields[$headers{'G'}]);
	my $N=safelog(0);
	push @matrixOutput,"$A\t$C\t$G\t$N\t$T\n";

	my $A=$fields[$headers{'A'}];
	my $T=$fields[$headers{'T'}];
	my $C=$fields[$headers{'C'}];
	my $G=$fields[$headers{'G'}];
	my $N=0;
	push @matrix,[$A,$C,$G,$N,$T];
      }
  }
my $matrixLength=@matrix;
my $cutoff=getCutoff();
my $alphabetSize=5;

@matrixOutput=();
for(my $i=0 ; $i<$matrixLength ; ++$i)
{
    my $row=$matrix[$i];
    my $nullRow=$nullWMM->[$i];
    my $n=@$row;
    for(my $j=0 ; $j<$n ; ++$j)
    {
	#print $row->[$j] ." -= ". $nullRow->[$j] ." (". $nullRow->[$j] .")\n";
	if($row->[$j]==0) 
	{
	    $row->[$j]="-inf";
	}
	else
	{
	    #$row->[$j]=safelog($row->[$j]/$nullRow->[$j]);
	    $row->[$j]=safelog($row->[$j]);
	    #print "\t\t" . $row->[$j] . "\n";
	}
    }
    push @matrixOutput,join("\t",@$row)."\n";
}

open(OUT,">$filestem.model") || die;
print OUT "WMM\nPROMOTER\n$cutoff\t$matrixLength\t$alphabetSize\n";
print OUT "$matrixLength\t$consensusOffset\t$consensusLength\n+\n";
my $matrix=join('',@matrixOutput);
print OUT "$matrix";
close(OUT);
close(IN);



sub safelog
  {
    my ($x)=@_;
    return (($x==0) ? "-inf" : int(100000*log($x)))/100000;
  }



sub buildGenerator
{
    my ($wmm)=@_;
    my @generator;
    my $wmmLen=@$wmm;
    for(my $i=0 ; $i<$wmmLen ; ++$i)
    {
	my $slot=[];
	push @generator,$slot;
	my $boundary=0;
	foreach my $elem (@{$wmm->[$i]})
	{
	    $boundary+=$elem;
	    push @$slot,$boundary;
	}
    }
    return \@generator;
}



sub getNullModel
{
    my ($wmm)=@_;
    my @nullModel;
    my $wmmLen=@$wmm;
    my $consensusEnd=$consensusOffset+$consensusLength;
    my $uniformSlot=[0.25,0.25,0.25,0,0.25]; # ACGNT
    
    for(my $i=0 ; $i<$wmmLen ; ++$i)
    {
	#if($i>=$consensusOffset && $i<$consensusEnd)
	#{
	#    push @nullModel,$wmm->[$i];
	#}
	#else
	{
	    my $realSlot=$wmm->[$i];
	    my $numZeros=0;
	    foreach my $entry (@$realSlot)
	    {
		if($entry==0) {++$numZeros}
	    }
	    my $numNonzero=5-$numZeros;
	    my $p=1/$numNonzero;
	    my $uniformSlot=[];
	    foreach my $entry (@$realSlot)
	    {
		if($entry==0) {push @$uniformSlot,0}
		else {push @$uniformSlot,$p}
	    }
	    push @nullModel,$uniformSlot;
	}
    }
    return \@nullModel;    
}



sub generate
{
    my ($model)=@_;
    my @strings;
    for(my $i=0 ; $i<$NUM_STRINGS ; ++$i)
    {
	my $string;
	for(my $j=0 ; $j<$matrixLength ; ++$j)
	{
	    my $slot=$model->[$j];
	    my $r=rand(1);
	    my $k;
	    for($k=0 ; $k<5 ; ++$k)
	    {
		my $boundary=$slot->[$k];
		if($r<$boundary) {last};
	    }
	    if($k==0) {$string.="A"}
	    elsif($k==1) {$string.="C"}
	    elsif($k==2) {$string.="G"}
	    elsif($k==3) {$string.="N"}
	    elsif($k==4) {$string.="T"}
	}
	push @strings,$string;
    }
    return \@strings;
}



sub scoreStrings
{
    my ($strings,$label,$scores)=@_;
    foreach my $string (@$strings)
    {
	my $score;
	for(my $i=0 ; $i<$matrixLength ; ++$i)
	{
	    my $slot=$matrix[$i];
	    my $base=substr($string,$i,1);
	    if($base eq "A") {$score+=safelog($slot->[0])}
	    elsif($base eq "C") {$score+=safelog($slot->[1])}
	    elsif($base eq "G") {$score+=safelog($slot->[2])}
	    elsif($base eq "N") {$score+=safelog($slot->[3])}
	    else {$score+=safelog($slot->[4])}
	}
	push @$scores,[$score,$label];
	#print "$label)$score : $string\n";
    }
}



sub getCutoff
{
    # FIRST, BUILD A STRING GENERATOR FROM THIS MODEL
    my $generator=buildGenerator(\@matrix);

    # BUILD A NULL MODEL USING BACKGROUND PROBABILITIES
    $nullWMM=getNullModel(\@matrix);

    foreach my $slot (@$nullWMM)
    {
	my $line=join("\t",@$slot);
	#print "$line\n";
    }

    # BUILD A STRING GENERATOR FROM THE NULL MODEL
    my $nullGenerator=buildGenerator($nullWMM);

    # GENERATE POSITIVE AND NEGATIVE STRINGS
    my $positiveStrings=generate($generator);
    my $negativeStrings=generate($nullGenerator);

    # SCORE ALL STRINGS WITH MODEL
    my @scores;
    scoreStrings($positiveStrings,1,\@scores);
    scoreStrings($negativeStrings,0,\@scores);

    # IDENTIFY CUTOFF FROM DESIRED SENSITIVITY
    open(SCORES,">$filestem.scores") || die;
    open(PRECRECALL,">$filestem.prec-recall") || die;
    @scores=sort {$a->[0] <=> $b->[0]} @scores;
    my ($pos,$neg);
    foreach my $elem (@scores) {if($elem->[1]){++$pos}else{++$neg}}
    my ($TP,$TN,$FP,$FN);
    my $numScores=@scores;
    my $selectedCutoff;
    for(my $i=1 ; $i<$numScores ; ++$i)
    {
	my $elem=$scores[$i];
	my $cutoff=$elem->[0];
	if($elem->[1]) {++$FN} else {++$TN}
	$TP=$pos-$FN;
	$FP=$neg-$TN;
	my $sens=$TP/$pos;
	my $spec=$TP/($TP+$FP);
	if((!defined($selectedCutoff) || $cutoff>$selectedCutoff) && 
	   $sens>=$minSensitivity)
	{
	    $selectedCutoff=$cutoff;
	}

	print PRECRECALL "$spec $sens\n";

	print SCORES "$elem->[1] $elem->[0]\n";
    }
    close(PRECRECALL);
    close(SCORES);
    return $selectedCutoff;
}

