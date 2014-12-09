#!/usr/bin/perl
use strict;
use FastaReader;
use ProgramName;

my $FALSE=0;
my $TRUE=1;
my $NUM_LEFT_SKIP_BUCKETS=1;
my $SMOOTH_LEFT_MARGIN=$FALSE;

my $name=ProgramName::get();
my $usage="
$name <*.fasta> <bin-size> <smoothing-window-size> 
          <smoothing-iterations> <#maxima> <maxima-multiplier>
";
die "$usage\n" unless @ARGV==6;
my ($filename,$binSize,$windowSize,$numIterations,$numMaxima,
    $maxFactor)=@ARGV;

$filename=~/(.*).fasta/ || die "$filename must have a .fasta extension\n";
my $outfile="$1.distr";

# Read fasta file & collect lengths
my (%counts,$total,@stack,$max);
my $reader=new FastaReader($filename);
while(1)
  {
    my ($defline,$sequence)=$reader->nextSequence();
    last unless defined $defline;
    my $length=length($sequence);
    push @stack,$length;
    if(!defined($max) || $length>$max) {$max=$length}
  }
my $n=@stack;

# Build histogram
my (@histogram,$largestBin,$total);
foreach my $x (@stack)
  {
    my $bin=int($x/$binSize);
    ++@histogram[$bin];
    if($bin>$largestBin){$largestBin=$bin}
    ++$total;
  }

# Normalize histogram & use set of pairs
for(my $i=0 ; $i<=$largestBin ; ++$i) 
  {
    my $count=$histogram[$i];
    my $p=$count/$total;
    $histogram[$i]=[$i*$binSize,$p];
  }

# Identify maxima & freeze them so they keep their original height
# even after smoothing
my @sorted=sort {$b->[1] <=> $a->[1]} @histogram;
my %maxima;
for(my $i=0 ; $i<$numMaxima ; ++$i)
  {
    $maxima{$sorted[$i]->[0]}=1;
    print STDERR "setting ". $sorted[$i]->[0] . " as a max\n";
  }

# Smooth histogram
for(my $i=0 ; $i<$numIterations ; ++$i)
  {
    @histogram=@{smooth($windowSize,\@histogram)};
  }

# Compute smoothness
my (@diffs,$prevDiff,$reversals);
for(my $i=1 ; $i<$n ; ++$i)
  {
    next unless $histogram[$i] && $histogram[$i-1];
    my ($thisX,$thisValue)=@{$histogram[$i]};
    my ($prevX,$prevValue)=@{$histogram[$i-1]};
    my $diff=$thisValue-$prevValue;
    #print "$prevValue $thisValue $diff\n";
    push @diffs,$diff;
    if(defined($prevDiff))
      {
	my $prevSign=sign($prevDiff);
	my $sign=sign($diff);
	if($prevSign!=$sign) {++$reversals}
      }
    $prevDiff=$diff;
  }
my $smoothness=($reversals ? 1/$reversals : 1);
print STDERR "SMOOTHNESS: $smoothness\n";

# Output histogram
for(my $i=0 ; $i<=$largestBin ; ++$i)
  {
    my $pair=$histogram[$i];
    next unless defined $pair;
    my ($x,$y)=@$pair;
    print "$x $y\n";
  }




#-------------------------------------------------
sub sign
  {
    my ($x)=@_;
    return ($x<0 ? -1 : 1);
  }
sub smooth
  {
    my ($windowSize,$histogram)=@_;
    my $n=@$histogram;
    #--$windowSize;
    my $halfWindow=int($windowSize/2);
    my $otherHalf=$windowSize-$halfWindow; # in case it's an odd number
    #my $first=$windowSize*$NUM_LEFT_SKIP_BUCKETS;
    my $first=$halfWindow;
    my $last=$n-1-$otherHalf;
    my $newHistogram=[];

    # Handle the leftmost bins (too close to edge to use a full window)
    my $boundarySum;
    if($SMOOTH_LEFT_MARGIN)
      {
	for(my $i=0 ; $i<$windowSize ; ++$i) 
	  {
	    my $pair=$histogram->[$i];
	    my $y=(defined($pair) ? $pair->[1] : 0);
	    $boundarySum+=$y;
	  }
	my $boundaryAve=$boundarySum/$windowSize;
	for(my $i=0 ; $i<$first ; ++$i) 
	  {
	    $newHistogram->[$i]=$histogram->[$i];
	    $newHistogram->[$i]->[1]=$boundaryAve;
	  }
      }
    else
      {
	for(my $i=0 ; $i<$first ; ++$i) 
	  {
	    $newHistogram->[$i]=$histogram->[$i];
	  }
      }

    # Handle the rightmost bins (too close to edge to use a full window)
    $boundarySum=0;
    for(my $i=$last+1 ; $i<$n ; ++$i) 
      {
	my $pair=$histogram->[$i];
	my $y=(defined($pair) ? $pair->[1] : 0);
	$boundarySum+=$y;
      }
    my $boundaryAve=$boundarySum/$windowSize;
    for(my $i=$last+1 ; $i<$n ; ++$i) 
      {
	$newHistogram->[$i]=$histogram->[$i];
	$newHistogram->[$i]->[1]=$boundaryAve;
      }

    # Handle the main part of the histogram
    for(my $i=$first ; $i<=$last ; ++$i)
      {
	my $pair=$histogram->[$i];
	my ($x,$y)=@$pair;
	if(!$maxima{$x})
	  {
	    for(my $j=0 ; $j<$halfWindow ; ++$j)
	      {
		my $pair=$histogram->[$i-1-$j];
		my ($leftX,$leftY)=(defined($pair) ? @$pair : (0,0));
		$y+=$leftY;
	      }
	    for(my $j=0 ; $j<$otherHalf ; ++$j)
	      {
		my $pair=$histogram->[$i+1+$j];
		my ($rightX,$rightY)=(defined($pair) ? @$pair : (0,0));
		$y+=$rightY;
	      }
	    $y/=($windowSize+1);
	  }
	else 
	  {
	    $y*=$maxFactor;
	    print STDERR "BOOSTING $x by $maxFactor\n";
	  }
	$newHistogram->[$i]=[$x,$y];
      }

    # Re-normalize
    my $sum=0;
    for(my $i=0 ; $i<$n ; ++$i) {$sum+=$newHistogram->[$i]->[1]}
    for(my $i=0 ; $i<$n ; ++$i) {$newHistogram->[$i]->[1]/=$sum}
    return $newHistogram;
  }
