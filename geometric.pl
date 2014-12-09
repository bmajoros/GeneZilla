#!/usr/bin/perl
use strict;

die "$0 <mean>\n" unless @ARGV==1;
my $mean=shift @ARGV;

my $q=1/$mean;

my $sum=0;
for(my $i=1 ; $i<10000 ; ++$i)
  {
    my $p=P($i);
    print "$i $p\n";
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
