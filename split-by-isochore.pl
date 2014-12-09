#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use FastaReader;
use FileHandle;
use FastaWriter;
use ProgramName;

my $MARGIN_AROUND_GENES=1000;

# PROCESS THE COMMAND LINE
my $name=ProgramName::get();
my $usage="$name <*.gff> <*.fasta> <out-dir> <#isochores> 0 10 20 ..etc.. 80 100";
die "$usage\n" unless @ARGV>=6;
my $inGff=shift @ARGV;
my $inFasta=shift @ARGV;
my $outDir=shift @ARGV;
my $numIsochores=shift @ARGV;
my $zero=shift @ARGV;
die "First isochore boundary must be 0" unless $zero==0;
my @boundaries;
for(my $i=0 ; $i<$numIsochores ; ++$i)
  {
    my $boundary=shift @ARGV;
    push @boundaries,$boundary;
  }
die "Last boundary must be 100" unless $boundaries[@boundaries-1]==100;

# LOAD THE GFF FILE & BIN BY SUBSTRATE
my $gffReader=new GffTranscriptReader();
my $transcripts=$gffReader->loadGFF($inGff);
my %transcriptsBySubstrate;
my $numTranscripts=@$transcripts;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    push @{$transcriptsBySubstrate{$transcript->getSubstrate()}},$transcript;
  }

# CREATE OUTPUT FILES
my $lowerBound=0;
my (@fastaFiles,@gffFiles);
for(my $i=0 ; $i<$numIsochores ; ++$i)
  {
    my $upperBound=$boundaries[$i];
    my $outFastaName="$outDir/iso$lowerBound-$upperBound.fasta";
    my $outGffName="$outDir/iso$lowerBound-$upperBound.gff";
    my $fastaFile=new FileHandle(">$outFastaName");
    my $gffFile=new FileHandle(">$outGffName");
    push @fastaFiles,$fastaFile;
    push @gffFiles,$gffFile;
    $lowerBound=$upperBound;
  }
my $fastaWriter=new FastaWriter;
open(HIST,">gc.txt") || die;

# PROCESS EACH SEQUENCE IN THE MULTI-FASTA FILE
my @binCounts;
my $fastaReader=new FastaReader($inFasta);
my $newSubstrateId=1000;
while(1)
  {
    # READ NEXT SEQUENCE & PARSE OUT SUBSTRATE ID
    my ($defline,$seqRef)=$fastaReader->nextSequenceRef();
    last unless defined $defline;
    $defline=~/^\s*>\s*(\S+)/ || die "can't parse $defline";
    my $substrateId=$1;
    my $seqLen=length $$seqRef;

    # PROCESS EACH TRANSCRIPT ON THIS SUBSTRATE
    my $localTranscripts=$transcriptsBySubstrate{$substrateId};
    next unless defined $localTranscripts;
    my $numLocal=@$localTranscripts;
    for(my $i=0 ; $i<$numLocal ; ++$i)
      {
	my $transcript=$localTranscripts->[$i];
	my $transcriptId=$transcript->getID();
	my $begin=$transcript->getBegin();
	my $end=$transcript->getEnd();
	$begin=($begin<1000 ? 0 : $begin-1000);
	$end=($end+1000<$seqLen ? $end+1000 : $seqLen);
	my $subseq=substr($$seqRef,$begin,$end-$begin);
	my $GC=getGCcontent($subseq,$begin,$end);
	if($GC<0) {die "Bad sequence ($begin,$end) for transcript $transcriptId : $subseq"}
	print HIST "$GC\n";
	$transcript->shiftCoords(-$begin);### 2/27/04
	$transcript->setSubstrate($newSubstrateId);
	my $gff=$transcript->toGff();
	for(my $i=0 ; $i<$numIsochores ; ++$i)
	  {
	    if($GC<$boundaries[$i])
	      {
		my $fastaFile=$fastaFiles[$i];
		my $gffFile=$gffFiles[$i];
		++$binCounts[$i];
		$fastaWriter->addToFasta(">$newSubstrateId",
					 $subseq,$fastaFile);
		++$newSubstrateId;
		print $gffFile $gff;
		last;
	      }
	  }
      }
  }
open(BINCOUNT,">bincounts.txt") || die;
print BINCOUNT "G+C%   \t#sequences\n";
my $lowerBound=0;
for(my $i=0 ; $i<$numIsochores ; ++$i)
  {
    my $upperBound=$boundaries[$i];
    print BINCOUNT "$lowerBound\%\-$upperBound\%  \t$binCounts[$i]\n";
    $lowerBound=$upperBound;
  }
close(BINCOUNT);

# CLEAN UP
close(HIST);
for(my $i=0 ; $i<$numIsochores ; ++$i)
  {
    my $fastaFile=$fastaFiles[$i];
    my $gffFile=$gffFiles[$i];
    close($fastaFile);
    close($gffFile);
  }

#----------------------------------------------------------------
sub getGCcontent
  {
    my ($subseq,$begin,$end)=@_;
    my $GC=($subseq=~s/([GC])/$1/g);
    my $ATGC=($subseq=~s/([ATGC])/$1/g);
    if($ATGC==0) {return -1}
    $GC/=$ATGC;
    return 100*$GC;
  }


