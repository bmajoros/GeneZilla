#!/usr/bin/perl
use strict;
use GffReader;
use FastaReader;
use Translation;
use FastaWriter;

$0=~/([^\/]+)\s*$/;
die "
$1 <*.gff> <*.fasta> <donor-consensuses> <acceptor-consensuses> 
         <donors-out.fasta> <acceptors-out.fasta>

  example:
     $1 splice.gff contigs.fasta GT,GC AG donors.fasta acceptors.fasta

"
 unless @ARGV==6;
my ($gffFile,$inFasta,$donorCons,$acceptorCons,$donorsFile,
    $acceptorsFile)=@ARGV;

my @donorCons=split/[,\.-]/,$donorCons;
my @acceptorCons=split/[,\.-]/,$acceptorCons;
my (%donorCons,%acceptorCons);
foreach my $d (@donorCons) {$donorCons{$d}=1}
foreach my $a (@acceptorCons) {$acceptorCons{$a}=1}

my $gffReader=new GffReader;
my $features=$gffReader->loadGFF($gffFile);
my $n=@$features;
my %bySubstrate;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $feature=$features->[$i];
    my $substrateId=$feature->getSubstrate();
    push @{$bySubstrate{$substrateId}},$feature;
  }

open(DONORS,">$donorsFile") || die "can't write to file: $donorsFile\n";
open(ACCEPTORS,">$acceptorsFile") 
  || die "can't write to file: $acceptorsFile\n";
my $fastaWriter=new FastaWriter;

my $fastaReader=new FastaReader($inFasta);
while(1)
  {
    my ($def,$seqRef)=$fastaReader->nextSequenceRef();
    last unless defined $def;
    $def=~/^\s*>\s*(\S+)/ || die "can't parse defline: $def\n";
    my $substrateId=$1;
    my $substrateLength=length $$seqRef;
    my $features=$bySubstrate{$substrateId};
    next unless $features;
    my $n=@$features;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $feature=$features->[$i];
	my $begin=$feature->getBegin()-80;
	my $end=$feature->getEnd()+80;
	my $strand=$feature->getStrand();
	my $type=$feature->getType();
	if($begin<0 || $end>$substrateLength) {next}
	my $feature=substr($$seqRef,$begin,$end-$begin);
	if($strand eq "-") 
	  {$feature=Translation::reverseComplement(\$feature)}
	my $consensus=substr($feature,80,2);
	my $hash=$type eq "donor" ? \%donorCons : \%acceptorCons;
	if(!$hash->{$consensus}) 
	  {print "WARNING! $consensus is not a valid $type consensus.\n";
	   next}
	my $defline=">$type /begin=$begin /end=$end /strand=$strand";
	my $handle=$type eq "donor" ? \*DONORS : \*ACCEPTORS;
	$fastaWriter->addToFasta($defline,$feature,$handle);
      }
  }

close(DONORS);
close(ACCEPTORS);
