#!/usr/bin/perl
use strict;

my $usage="$0 <transfac-matrix> <consensus-offset> <consensus-length>
   <#strings> <pos-out.fasta> <neg-out.fasta>";
die "$usage\n" unless @ARGV==6;
my ($filename,$consensusOffset,$consensusLength,$numStrings,$posOut,
    $negOut)=@ARGV;

# ### HACK -- USE HUMAN NONCODING PROBABILITIES FOR NULL MODEL
my $pA=exp(-1.37373);
my $pC=exp(-1.45735);
my $pG=exp(-1.43431);
my $pT=exp(-1.28926);
my $pSum=$pA+$pC+$pG+$pT;
$pA/=$pSum;
$pC/=$pSum;
$pT/=$pSum;
$pG/=$pSum;

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
my $generator=buildGenerator(\@matrix);
open(OUT,">$posOut") || die;
generate($generator,\*OUT);
close(OUT);

my $nullModel=getNullModel(\@matrix);
my $nullGenerator=buildGenerator($nullModel);
open(OUT,">$negOut") || die;
generate($nullGenerator,\*OUT);
close(OUT);


#------------------------------------------------------------------------
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



sub generate
{
    my ($model,$handle)=@_;
    my @strings;
    for(my $i=0 ; $i<$numStrings ; ++$i)
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
	    else {$string.="T"}
	}
	print $handle ">$i (generated from TRANSFAC matrix)\n$string\n";
    }
}



sub getNullModel
{
    my ($wmm)=@_;
    my @nullModel;
    my $wmmLen=@$wmm;
    my $consensusEnd=$consensusOffset+$consensusLength;
    #my $uniformSlot=[0.25,0.25,0.25,0,0.25]; # ACGNT
    my $uniformSlot=[$pA,$pC,$pG,0,$pT];
    
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
