#!/usr/bin/perl
use strict;
use ProgramName;
use TempFilename;
$|=1;

# Process the command line
my $name=ProgramName::get();
my $usage="$name <*.gff> <chrom-dir> <margin-size>";
die "$usage\n" unless @ARGV==3;
my ($gffFile,$chromDir,$marginSize)=@ARGV;
if(-e "chunks") {die "Please rm chunks directory first\n"}
if(-e "iso0-100.fasta") {die "Please rm iso0-100.fasta\n"}
if(-e "iso0-100.gff") {die "Please rm iso0-100.gff\n"}

# Get the list of chromosomes
my $ls=`ls $chromDir/*.cof`;
my @ls=split/\s+/,$ls;
my $numChroms=@ls;
die "Chromosomes must be in *.cof format! (no defline or whitespace)\n" 
  unless $numChroms>0;

# Process each chromosome separately
my $nextChunkId=1;
foreach my $chromFile (@ls)
  {
    if(-e "chunks") {system("rm -r chunks")}

    # Use "make-test-set.pl" to extract the region around each transcript
    $chromFile=~/([^\/]+)\.cof/;
    my $substrateId=$1;
    print STDERR "Processing $substrateId...\n";
    my $cmd="make-test-set.pl -c $substrateId $gffFile $chromFile 1000000 $marginSize";
    print STDERR "$cmd\n";
    system($cmd);

    # Need to change the deflines and GFF substrate IDs
    my $ls=`ls chunks/*.fasta`;
    my @ls=split/\s+/,$ls;
    foreach my $chunkFile (@ls)
      {
	$chunkFile=~/([^\/]+)\.fasta/;
	my $chunkId=$1;
	my $gffFile="chunks/$chunkId.gff";
	my $cmd="fasta-change-defline.pl $chunkFile $nextChunkId";
	system($cmd);
	$cmd="gff-change-substrate.pl $gffFile $nextChunkId";
	system($cmd);
	++$nextChunkId;
      }

    # Gather the chunks from this chromosome into a single file
    $cmd="cat chunks/*.fasta > tmp.1";
    system($cmd);
    if(-e "iso0-100.fasta") 
      {system("cat tmp.1 iso0-100.fasta > tmp.2 ; mv tmp.2 tmp.1")}
    system("mv tmp.1 iso0-100.fasta");

    $cmd="cat chunks/*.gff > tmp.1";
    system($cmd);
    if(-e "iso0-100.gff") 
      {system("cat tmp.1 iso0-100.gff > tmp.2 ; mv tmp.2 tmp.1")}
    system("mv tmp.1 iso0-100.gff");
  }








