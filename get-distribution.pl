#!/usr/bin/perl
use strict;
use FastaReader;





die "OBSOLETE!  Use get-duration-distr.pl!\n";






my $usage="$0 <*.fasta> <smoothing-window-size-or-0> <smoothing-iterations>";
die "$usage\n" unless @ARGV==3;
my ($filename,$windowSize,$iterations)=@ARGV;

$filename=~/(.*).fasta/ || die "$filename must have a .fasta extension\n";
my $outfile="$1.distr";

my (%counts,$total,@stack,$max);
my $reader=new FastaReader($filename);
while(1)
  {
    my ($defline,$sequence)=$reader->nextSequence();
    last unless defined $defline;
    #next unless rand(1)<0.2;
    my $length=length($sequence);
    push @stack,$length;
    if(!defined($max) || $length>$max) {$max=$length}
  }
my $n=@stack;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $length=$stack[$i];
    ++$counts{int($length)};
    ++$total;
  }

my @histogram;
my @lengths=keys %counts;
@lengths=sort {$a<=>$b} @lengths;
my $numLengths=@lengths;
for(my $i=0 ; $i<$numLengths ; ++$i)
  {
    my $length=$lengths[$i];
    my $count=$counts{$length};
    my $p=$count/$total;
    push @histogram,[$length,$p];
  }

open(OUT,">$outfile") || die "can't create $outfile\n";
if($windowSize>0 && $iterations>0)
  {
    my $histogram=smooth($windowSize,\@histogram);
    for(my $i=1 ; $i<$iterations ; ++$i)
      {$histogram=smooth($windowSize,$histogram)}
    printHistogram($histogram);
  }
else
  {
    printHistogram(\@histogram);
  }
close(OUT);
print STDERR "sample size $numLengths\n";

#-------------------------------------------------
sub printHistogram
  {
    my ($histogram)=@_;
    my $n=@$histogram;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $pair=$histogram->[$i];
	my ($x,$y)=@$pair;
	print OUT "$x $y\n";
      }
  }
#-------------------------------------------------
sub smooth
  {
    my ($windowSize,$histogram)=@_;
    my $n=@$histogram;
    --$windowSize;
    my $halfWindow=int($windowSize/2);
    my $otherHalf=$windowSize-$halfWindow; # in case it's an odd number
    my $first=$halfWindow;
    my $last=$n-1-$otherHalf;
    my $newHistogram=[];
    for(my $i=0 ; $i<$first ; ++$i) {$newHistogram->[$i]=$histogram->[$i]}
    for(my $i=$last+1 ; $i<$n ; ++$i) {$newHistogram->[$i]=$histogram->[$i]}
    for(my $i=$first ; $i<=$last ; ++$i)
      {
	my $pair=$histogram->[$i];
	my ($x,$y)=@$pair;
	for(my $j=0 ; $j<$halfWindow ; ++$j)
	  {
	    my $debug=$i-1-$j;
	    my ($leftX,$leftY)=@{$histogram->[$i-1-$j]};
	    $y+=$leftY;
	  }
	for(my $j=0 ; $j<$otherHalf ; ++$j)
	  {
	    my $pair=$histogram->[$i+1+$j];
	    my ($rightX,$rightY)=(defined($pair) ? @$pair : (0,0));
	    $y+=$rightY;
	  }
	$y/=($windowSize+1);
	$newHistogram->[$i]=[$x,$y];
      }
    return $newHistogram;
  }



