#!/usr/bin/perl
use strict;
use Progress;

my $usage="$0 <*.iso>";
die "$usage\n" unless @ARGV==1;
my ($isoFile)=@ARGV;

if(!-e "out") {`mkdir out`}
my $ls=`ls chunks/*.fasta`;
my @files=split/\s+/,$ls;
my $numChunks=@files;

system("rm out/*.gff");

my $progress=new Progress;
$progress->start($numChunks);
foreach my $i (1..$numChunks)
  {
    my $chunk="chunks/$i.fasta";
    next unless -e $chunk;
    my $command="genezilla/genezilla $isoFile $chunk > out/$i.gff";
    print "$command\n";
    system($command);
    my ($timeLeft,$percentDone)=$progress->getProgress($i);

    my $msg="$percentDone\% done.  $timeLeft remaining.";
    my $msgLen=length $msg;
    my $bar='='x$msgLen;
    print "\n$bar\n$msg\n$bar\n\n";
  }
