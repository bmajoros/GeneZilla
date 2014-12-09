#!/usr/bin/perl
use strict;
use SummaryStats;

my (%signals,%contents,$totalLength);
while(<STDIN>)
  {
    if(/signal\s+(\S+)\s+(\S+)/)
      {
	my ($signalType,$score)=($1,$2);
	$signalType=~s/-//g;
	push @{$signals{$signalType}},$score;
      }
    elsif(/content\s+(\S+)\s+(\S+)\s+(\S+)/)
      {
	my ($contentType,$score,$length)=($1,$2,$3);
	$contentType=~s/NEG-//g;
	if(!defined($contents{$contentType})) {$contents{$contentType}={}}
	my $record=$contents{$contentType};
	push @{$record->{scores}},$score/$length;
	push @{$record->{lengths}},$length;
      }
    elsif(/TOTAL LENGTH: (\d+)/) {$totalLength=$1}
  }

my @signalTypes=keys %signals;
my @contentTypes=keys %contents;

foreach my $signalType (@signalTypes)
  {
    my $record=$signals{$signalType};
    my ($mean,$sd,$min,$max)=SummaryStats::summaryStats($record);
    $mean=int(10*$mean+0.5)/10;
    $sd=int(10*$sd+0.5)/10;
    my $n=@$record;
    my $density=$n/$totalLength;
    print "$signalType: $mean +/- $sd (N=$n; density=$density)\n"
  }

foreach my $contentType (@contentTypes)
  {
    my $record=$contents{$contentType};
    my ($meanScore,$sdScore,$minScore,$maxScore)=
      SummaryStats::summaryStats($record->{scores});
    my $n=@{$record->{scores}};
    $meanScore=int(10*$meanScore+0.5)/10;
    my ($meanLen,$sdLen,$minLen,$maxLen)=
      SummaryStats::summaryStats($record->{lengths});
    $meanLen=int(10*$meanLen+0.5)/10;
    print "$contentType SCORES: $meanScore (N=$n)\n";
    print "\t\tLENGTHS: $meanLen (N=$n)\n"
  }

