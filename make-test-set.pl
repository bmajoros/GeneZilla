#!/usr/bin/perl
use strict;
use FastaReader;
use CompiledFasta;
use FastaWriter;
use GffTranscriptReader;
use Translation;
use FileSize;
use ProgramName;
use Getopt::Std;
$|=1;

our $opt_c;
getopts('c:');

my $name=ProgramName::get();
my $usage="$name [-c substrate_id] <*.gff> <*.fasta> <max#> <margin-size> TAG,TGA,TAA

where -c indicates that the fasta file is compiled (*.cof)
         (compiled fasta files contain no defline and no whitespace)
";
die "$usage\n" unless @ARGV==5;
my ($gffFile,$fastaFile,$maxN,$margin,$stops)=@ARGV;
die "$gffFile does not exist\n" unless -e $gffFile;
die "$fastaFile does not exist\n" unless -e $fastaFile;

if(!-e "chunks") {system("mkdir chunks")}

my @stops=split/\s+/,$stops;
my %stops;
foreach my $stop (@stops) {$stops{$stop}=1}

##################### LOAD GFF #######################
my $gffReader=new GffTranscriptReader;
$gffReader->setStopCodons(\%stops);
my $transcripts=$gffReader->loadGFF($gffFile);
my $numTranscripts=@$transcripts;

# Index by contig ID and add padding around each gene
my %records;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    my $begin=$transcript->getBegin();
    my $end=$transcript->getEnd();
    die unless $begin<$end;
    next unless $begin>1;
    $begin-=$margin;
    $end+=$margin;
    my $substrate=$transcript->getSubstrate();
    if(defined($opt_c)) {next unless $substrate eq $opt_c}
    my $strand=$transcript->getStrand();
    push @{$records{"\L$substrate"}},[$begin,$end-$begin,$strand,$transcript];
  }

# Eliminate any gene having alternative splicing
my @substrates=keys %records;
my $n=@substrates;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $substrate=$substrates[$i];
    my $recordsOnSubstrate=$records{$substrate};
    next unless defined $recordsOnSubstrate;
    my $nn=@$recordsOnSubstrate;
    for(my $j=0 ; $j<$nn ; ++$j)
      {
	my $record1=$recordsOnSubstrate->[$j];
	my $transcript1=$record1->[3];
	my $begin1=$transcript1->getBegin();
	my $end1=$transcript1->getEnd();

	for(my $k=0 ; $k<$nn ; ++$k)
	  {
	    next if $k==$j;
	    my $record2=$recordsOnSubstrate->[$k];
	    my $transcript2=$record2->[3];
	    my $begin2=$transcript2->getBegin();
	    my $end2=$transcript2->getEnd();
	    if($begin2<$end1 && $end2>$begin1)
	      {$transcript2->{alt}=1}
	  }
      }
  }

################ PROCESS FASTA FILE #################
my $fastaReader;
if($opt_c)
  {$fastaReader=new CompiledFasta($fastaFile)}
else
  {$fastaReader=new FastaReader($fastaFile)}

my $fastaWriter=new FastaWriter;
my $nextId=1;
my %enumeration;
my %debugHash;
my $count;
while(1)
  {
    last unless $count<$maxN;
    my ($defline,$seq);
    if($opt_c)
      {$defline=">$opt_c"}
    else
      {($defline,$seq)=$fastaReader->nextSequenceRef()}
    last unless defined $defline;
    my $seqLen=($opt_c ? FileSize::fileSize($fastaFile) : length $$seq);
    $defline=~/>(\S+)/ || die $defline;
    my $bac="\L$1";
    my $records=$records{$bac};
    if(!defined($records))
       {
	 if($opt_c) {last}
	 next;
       }
    my $n=@$records;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	last unless $count<$maxN;
	my $record=$records->[$i];
	my ($begin,$len,$strand,$transcript)=@$record;
	next if($transcript->{alt});

	if($begin<0) {$begin=0}
	my $end=$begin+$len;#+3;
	if($end>$seqLen) {$end=$seqLen}
	$len=$end-$begin;

	my $transcriptId=$transcript->getID();
	if($transcriptId=~/(.*);/) {$transcriptId=$1}
	if(!defined($enumeration{$transcriptId})) 
	  {$enumeration{$transcriptId}=$nextId++}

	$debugHash{$enumeration{$transcriptId}}=$transcriptId;
	my $d=$enumeration{$transcriptId};

	$transcriptId=$enumeration{$transcriptId};
	my $outfasta="$transcriptId.fasta";
	my $outgff="$transcriptId.gff";
	open(OUTFASTA,">chunks/$outfasta") || die "can't create $outfasta";
	open(OUTGFF,">chunks/$outgff") || die "can't create $outgff";

	my $string;
	if($opt_c)
	  {$string=$fastaReader->load($begin,$len)}
	else
	  {$string=substr($$seq,$begin,$len)}
	my $substrate=$transcript->getSubstrate();
	$fastaWriter->addToFasta(">$transcriptId ($substrate)",$string,
				 \*OUTFASTA);
	$transcript->setSubstrate($transcriptId);

	# Shift coords and output GFF
	my $numExons=$transcript->numExons();
	my $delta=-$begin;
	for(my $j=0 ; $j<$numExons ; ++$j)
	  {
	    my $exon=$transcript->getIthExon($j);
	    $exon->shiftCoords($delta);
	    my $gff=$exon->toGff();
	    print OUTGFF "$gff";
	  }
	close(OUTGFF);
	close(OUTFASTA);
	++$count;
      }
    if($opt_c) {last}
  }

