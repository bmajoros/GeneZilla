#!/usr/bin/perl
##########################################################################
# This program estimates transition probabilities for the GHMM
# bmajoros@tigr.org
##########################################################################
use strict;
use GffTranscriptReader;

my $usage="$0 <*.gff>";
die "$usage\n" unless @ARGV==1;
my ($chunksGff)=@ARGV;

my $gffReader=new GffTranscriptReader();
my $transcripts=$gffReader->loadGFF($chunksGff);#"chunks.gff");

my ($AG_GT,$AG_TAG,$ATG_GT,$ATG_TAG);
my (@intronPhases,@GT,@AG,@initGT,@finalAG);
my $numTranscripts=@$transcripts;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    my $strand=$transcript->getStrand();
    my $numExons=$transcript->numExons();

    if($numExons==1) {++$ATG_TAG}
    else
      {
	++$ATG_GT;
	++$AG_TAG;
	my $numAG=$numExons-1;
	$AG_GT+=($numAG-1);

	# Now handle intron phases
	for(my $i=0 ; $i<$numExons ; ++$i)
	  {
	    my $exon=$transcript->getIthExon($i);
	    my $exonFrame=$exon->getFrame();
	    my $exonLength=$exon->getLength();
	    my $beginPhase=$exonFrame;
	    my $endPhase=($beginPhase+$exonLength)%3;
	    my $hasDonor=($i+1<$numExons);
	    my $hasAcceptor=($i>0);
	    if($hasDonor) {++$GT[$endPhase]}
	    if($hasAcceptor) {++$AG[$beginPhase]}
	    if($hasDonor && $hasAcceptor) 
	      {++$intronPhases[$beginPhase]->[$endPhase]}
	    if(!$hasAcceptor && $hasDonor) {++$initGT[$endPhase]} # initial exon
	    if($hasAcceptor && !$hasDonor) {++$finalAG[$beginPhase]} # final exon
	  }
      }
  }
my $AG=$AG_GT+$AG_TAG;
my $GT=$AG_GT+$ATG_GT;
my $ATG=$ATG_GT+$ATG_TAG;
my $TAG=$AG_TAG+$ATG_TAG;

my $P_ATG_GT=$ATG_GT/$ATG;
my $P_ATG_TAG=$ATG_TAG/$ATG;
my $P_AG_GT=$AG_GT/$AG;
my $P_AG_TAG=$AG_TAG/$AG;

my $P_mGT_mATG=$ATG_GT/$GT;
my $P_mGT_mAG=$AG_GT/$GT;
my $P_mTAG_mATG=$ATG_TAG/$TAG;
my $P_mTAG_mAG=$AG_TAG/$TAG;

my $numSignalTypes=24; # |{ATG,GT,AG,TAG,POLYA,PROMOTER,GT#,AG#}|=12, *2=24

print "$numSignalTypes\n";

# FORWARD STRAND:
print "PROMOTER -> ATG : 1\n";
print "ATG -> GT : $P_ATG_GT\n";
print "ATG -> TAG : $P_ATG_TAG\n";
print "GT -> AG : 1\n";
print "AG -> GT : $P_AG_GT\n";
print "AG -> TAG : $P_AG_TAG\n";
print "TAG -> POLYA : 0.5\n";
print "TAG -> ATG : 0.125\n";
print "TAG -> PROMOTER : 0.125\n";
print "TAG -> -TAG : 0.125\n";
print "TAG -> -POLYA : 0.125\n";
print "POLYA -> ATG : 0.25\n";
print "POLYA -> PROMOTER : 0.25\n";
print "POLYA -> -TAG : 0.25\n";
print "POLYA -> -POLYA : 0.25\n";
for(my $i=0 ; $i<3 ; ++$i)
  {
    print "GT$i -> AG$i : 1\n";
    my $AGi=$AG[$i];
    for(my $j=0 ; $j<3 ; ++$j)
      {
	my $AGi_GTj=$intronPhases[$i]->[$j];
	my $p=$AGi_GTj/$AGi;
	print "AG$i -> GT$j : $p\n";
      }
    my $GTi=$initGT[$i];
    my $ATG_GTi=$GTi/$ATG;
    print "ATG -> GT$i : $ATG_GTi\n";
    my $final_AGi=$finalAG[$i];
    my $AGi_TAG=$final_AGi/$AGi;
    print "AG$i -> TAG : $AGi_TAG\n";
  }


# REVERSE STRAND:
print "-PROMOTER -> ATG : 0.25\n";
print "-PROMOTER -> PROMOTER : 0.25\n";
print "-PROMOTER -> -TAG : 0.25\n";
print "-PROMOTER -> -POLYA : 0.25\n";
print "-ATG -> -PROMOTER : 0.5\n";
print "-ATG -> ATG : 0.125\n";
print "-ATG -> PROMOTER : 0.125\n";
print "-ATG -> -TAG : 0.125\n";
print "-ATG -> -POLYA : 0.125\n";
print "-GT -> -ATG : $P_mGT_mATG\n";
print "-GT -> -AG : $P_mGT_mAG\n";
print "-AG -> -GT : 1\n";
print "-TAG -> -ATG : $P_mTAG_mATG\n";
print "-TAG -> -AG : $P_mTAG_mAG\n";
print "-POLYA -> -TAG : 1\n";
for(my $i=0 ; $i<3 ; ++$i)
  {
    print "-AG$i -> -GT$i : 1\n";
    my $GTi=$GT[$i];
    for(my $j=0 ; $j<3 ; ++$j)
      {
	my $AGj_GTi=$intronPhases[$j]->[$i];
	my $p=$AGj_GTi/$GTi;
	print "-GT$i -> -AG$j : $p\n";
      }
    my $init_GTi=$initGT[$i];
    my $GTi_ATG=$init_GTi/$GTi;
    print "-GT$i -> -ATG : $GTi_ATG\n";
    my $final_AGi=$finalAG[$i];
    my $TAG_AGi=$final_AGi/$TAG;
    print "-TAG -> -AG$i : $TAG_AGi\n";
  }
