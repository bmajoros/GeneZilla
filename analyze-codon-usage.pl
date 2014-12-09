#!/usr/bin/perl
use strict;
use FastaReader;
use Translation;
use Feature;
$|=1;

# PROCESS COMMAND LINE
die "$0 <*.fasta> <ORFs.gff> <min-ORF-length>\n" 
  unless @ARGV==3;
my ($fastaFile,$gffFile,$minLength)=@ARGV;

# LOAD ORF COORDINATES
print STDERR "Reading GFF file...\n";
my $features=loadGFF($gffFile,$minLength);

# PUT ORFS IN HASH TABLE BY SUBSTRATE ID
print STDERR "Indexing features...\n";
my ($totalSeq,%hash,@profiles,%centroid,$globalTotal);
my $numFeatures=@$features;
for(my $i=0 ; $i<$numFeatures ; ++$i)
  {
    my $feature=$features->[$i];
    my $substrateId=$feature->getSubstrate();
    my $length=$feature->getLength();
    push @{$hash{$substrateId}},$feature;
  }

# PROCESS ALL SUBSTRATE SEQUENCES
print STDERR "Processing contigs file...\n";
my $fastaReader=new FastaReader($fastaFile);
while(1)
  {
    my ($def,$seq)=$fastaReader->nextSequence();
    last unless defined $def;
    $def=~/^>(\S+)/ || die "Can't parse defline: $def\n";
    my $substrateId=$1;
    my $features=$hash{$substrateId};
    next unless defined $features;
    my $length=length $seq;
    $totalSeq+=$length;
    print STDERR "TOTAL SEQUENCE PROCESSED: $totalSeq\n";
    my $numFeatures=@$features;
    for(my $i=0 ; $i<$numFeatures ; ++$i)
      {
	my $feature=$features->[$i];
	my $begin=$feature->getBegin();
	my $end=$feature->getEnd();
	my $length=$end-$begin;
	my $subseq=substr($seq,$begin,$length);
	next if($subseq=~/N/);
	my $strand=$feature->getStrand();
	if($strand eq '-') {$subseq=Translation::reverseComplement(\$subseq)}
	my $phase=$feature->getFrame();
	my $profile={};
	processCodons($subseq,$phase,$length,$feature,$profile);
	push @profiles,{vector=>$profile, seq=>$subseq, phase=>$phase};
      }

    undef $hash{$substrateId};
    my @ids=keys %hash;
    last if(@ids==0);
  }

# COMPUTE CENTROID
my @codons=keys %centroid;
foreach my $codon (@codons)
  {$centroid{$codon}/=$globalTotal}

# NOW PROCESS ALL CODON-USAGE PROFILES
print STDERR "Processing codon-usage profiles...\n";
my $numProfiles=@profiles;
for(my $i=0 ; $i<$numProfiles ; ++$i)
  {
    my $profile=$profiles[$i];
    my $vector=$profile->{vector};
    my $distance=relativeEntropy($vector,\%centroid);
    print "$distance\n";

    #if($distance>=0.4)
    #  {
	#my $seq=$profile->{seq};
	#print "$distance $seq\n";
      #}
  }

#----------------------------------------------------------
sub processCodons
  {
    my ($seq,$phase,$length,$feature,$hash)=@_;
    my $begin=(3-$phase)%3;
    my $total;
    for(my $pos=$begin ; $pos+2<$length ; $pos+=3)
      {
	my $codon=substr($seq,$pos,3);
	++$hash->{$codon};
	++$total;

	++$centroid{$codon};
	++$globalTotal;
      }
    my @codons=keys %$hash;
    foreach my $codon (@codons)
      {$hash->{$codon}/=$total}
  }
#----------------------------------------------------------
sub relativeEntropy
  {
    my ($profileI,$profileJ)=@_;
    my $relEnt=0;
    my @codons=keys %$profileI;
    foreach my $codon (@codons)
      {
	###print "codon: $codon ";
	my $p=$profileI->{$codon};
	my $q=$profileJ->{$codon};
	###print "p=$p q=$q\n";
	next unless $p>0 && $q>0;
	$relEnt+=($p*log($p/$q));
      }
    return $relEnt;
  }
#----------------------------------------------------------
sub loadGFF
  {
    my ($gffFilename,$minLength)=@_;
    my @features;
    if($gffFilename=~/\.gz$/)
      {open(GFF,"cat $gffFilename|gunzip|") || die $gffFilename}
    else
      {open(GFF,$gffFilename) || die $gffFilename}
    while(<GFF>)
      {
	next unless $_=~/\S+/;
	next if $_=~/^\s*\#/;
	my @fields=split/\s+/,$_;
	next unless @fields>7;
	my $length=$fields[4]-$fields[3];
	next unless $length>=$minLength;
	my @additionalFields=splice(@fields,9);
	my $feature=new Feature($fields[0],$fields[1],$fields[2],
				$fields[3],$fields[4],$fields[5],
				$fields[6],$fields[7],$fields[8],
				\@additionalFields);
	push @features,$feature;
	
	###
	#last if @features>4;
	###
      }
    close(GFF);
    @features=sort {$a->{begin} <=> $b->{begin}} @features;
    return \@features;
  }
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
