#!/usr/bin/perl
use strict;

die "$0 <out-filestem>\n" unless @ARGV==1;
my ($filestem)=@ARGV;

my $outfile="$filestem.tar";
my $files="*.model *.top *.iso *.cfg *.distr *.trans chunks/*.gff chunks/*.fasta";
my $command="tar cf $outfile $files ; gzip $outfile";

print "$command\n";
system($command);

$outfile.=".gz";
print "Archive $outfile created\n";

