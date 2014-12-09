#!/usr/bin/perl
use strict;

my $browser="/home/bmajoros/genomics/browser/make-figure";
my $base=".";
my $bacsPerPage=5;
my $maxPages=1000;

if(!-d "ps") {system("mkdir ps")}
my $page=1;
for(my $i=1 ; $i<4000 ; $i+=$bacsPerPage, ++$page)
  {
    my $command="$browser 2";
    my $found=0;
    for(my $bac=0 ; $bac<$bacsPerPage ; ++$bac)
      {
	my $id=$i+$bac;
	my $true="$base/chunks/$id.gff";
	my $ghmm="$base/out/$id.gff";
	next unless -e $ghmm;
	$command.=" $true correct\#$id $ghmm predict\#$id";
	#$command.=" $true cDNA\#$id $ghmm GeneZilla\#$id";
	++$found;
      }
    next unless $found;
    $command.=" ps/page$page.ps";
    system($command);
    last unless $page<$maxPages;
  }

