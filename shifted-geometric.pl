#!/usr/bin/perl
use strict;

die "$0 <mean> <delta-x> <bucket-size>\n" unless @ARGV==3;
my $mean=shift @ARGV;
my $deltaX=shift @ARGV;
my $bucketSize=shift @ARGV;

my $q=1/$mean;

my $sum=0;
for(my $i=1 ; $i<10000 ; ++$i)
  {
    my $p=$bucketSize*P($i);
    my $x=$i+$deltaX;
    print "$x $p\n";
    my $d=$i*$p;
    $sum+=$d;
  }
#print "$sum\n";

sub P
  {
    my ($len)=@_;
    my $p=((1-$q)**($len-1))*$q;
    return $p;
  }
