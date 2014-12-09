#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <pos-class.score> <neg-class.scores>\n" unless @ARGV==2;
my ($posFile,$negFile)=@ARGV;

my $posScores=load($posFile);
my $negScores=load($negFile);

my $N=@$posScores+@$negScores;
my $correct=greater1($posScores)+@$negScores-greater1($negScores);
my $acc=$correct/$N;
print "$acc\n";

sub greater1
  {
    my ($scores)=@_;
    my $n=@$scores;
    my $count=0;
    for(my $i=0 ; $i<$n ; ++$i) {
      if($scores->[$i]>1) { ++$count }
    }
    return $count;
  }

sub load
  {
    my ($filename)=@_;
    my $scores=[];
    open(IN,$filename) || die "can't open $filename\n";
    while(<IN>) {
      chomp;
      if(/^\s*(\S+)/) { push @$scores,0+$1 }
    }
    close(IN);
    return $scores;
  }


