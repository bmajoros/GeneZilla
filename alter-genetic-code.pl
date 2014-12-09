#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;

# PROCESS COMMAND LINE
die "$0 <infile> <outfile>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

# READ SEQUENCES
my $seqs=FastaReader::readAllAndKeepDefs($infile);
my @seqs=values %$seqs;

# DUPLICATE TRAINING SEQUENCES AND PERFORM SUBSTITUTIONS IN THE
# SECOND HALF OF THE DUPLICATED SET
@seqs=substitute("ACC","TAG",\@seqs);
@seqs=substitute("TCC","TAA",\@seqs);

# WRITE OUTPUT
my $writer=new FastaWriter;
open(OUT,">$outfile") || die "can't create file: $outfile\n";
my $n=@seqs;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my ($def,$seq)=@{$seqs[$i]};
    $writer->addToFasta($def,$seq,\*OUT);
  }
close(OUT);

#---------------------------------------------------------------
sub substitute
  {
    my ($codon1,$codon2,$seqs)=@_;
    my @newSeqs=@$seqs;
    my $n=@$seqs;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $seq=$seqs->[$i];
	my $newSeq=performSubstitutions($seq,$codon1,$codon2);
	push @newSeqs,$newSeq;
      }
    return @newSeqs;
  }
#---------------------------------------------------------------
sub performSubstitutions
  {
    my ($defseq,$codon1,$codon2)=@_;
    my ($def,$seq)=@$defseq;
    $def=~/\/frame=(\d)/ || die "can't parse defline: $def";
    my $frame=$1;

    # Perform the actual (in-frame) substitution:
    my $begin=(3-$frame)%3;
    my $len=length $seq;
    for(my $pos=$begin ; $pos<$len ; $pos+=3)
      {
	if(substr($seq,$pos,3) eq $codon1)
	  {
	    substr($seq,$pos,3)=$codon2;
	    #print "substituted $codon2 for $codon1 at $pos\n";
	  }
      }

    return [$def,$seq];
  }
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

