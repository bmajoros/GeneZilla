#!/usr/bin/perl
use strict;

print "ngram\tphase 0    phase 1    phase 2\n";
print "-----\t-------    -------    -------\n";

$_=<STDIN>; chop;
if($_ ne "3P") {die "Model type is $_; not 3P\n"}

my %hash;
my $strand=1;
for(my $i=0 ; $i<6 ; ++$i)
  {
    my $MC=<STDIN>; chop $MC;
    die "$MC not MC!\n" unless $MC eq "MC";
    <STDIN>; # content type
    $_=<STDIN>;
    $_=~/(\d+)\s+(\d+)/;
    my ($n,$phase)=($1,$2);
    my $numModels=0+<STDIN>;
    for(my $j=0 ; $j<$numModels ; ++$j)
      {
	my $numElements=0+<STDIN>;
	for(my $i=0 ; $i<$numElements ; ++$i)
	  {
	    my $ngram=<STDIN>; chop $ngram;
	    my $len=length $ngram;
	    my $score=0+<STDIN>;
	    $score=int(100*$score+0.5)/100;
	    if($strand>0)
	      {
		$hash{$len}->{$ngram}->{$phase}=$score;
	      }
	  }
      }
    $strand*=-1;
  }

my @lengths=sort {$a <=> $b} keys %hash;
foreach my $len (@lengths)
  {
    my $hash2=$hash{$len};
    my @ngrams=sort keys %$hash2;
    foreach my $ngram (@ngrams)
      {
	if($ngram=~/N/) {next}
	print "$ngram \t";
	my $score0=pad($hash2->{$ngram}->{0},10);
	my $score1=pad($hash2->{$ngram}->{1},10);
	my $score2=pad($hash2->{$ngram}->{2},10);
	print "$score0 $score1 $score2\n";
      }
  }

sub pad
  {
    my ($n,$l)=@_;
    my $len=length $n;
    $n.=' 'x($l-$len);
    return $n;
  }
