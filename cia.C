/****************************************************************
 cia.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "CIA.H"
#include "cia.H"
#include "EdgeFactory.H"
#include "BOOM/FastaReader.H"
#include "BOOM/Constants.H"
using namespace std;

#ifndef EXPLICIT_GRAPHS
#error Please edit the file genezilla.H, uncomment the definition of EXPLICIT_GRAPHS, issue a "make clean", and recompile this project
#endif

#ifdef FORCE_SPECIFIC_SIGNALS
#error Please edit the file genezilla.H, comment out the definition of FORCE_SPECIFIC_SIGNALS, issue a "make clean", and recompile this project
#endif

static const char *PROGRAM_NAME="CIA";
static const char *VERSION="0.1";
Alphabet alphabet;
int frame; // ### CAUTION: this is required by older code; to be removed

void AppMain(int argc,char *argv[]);

int main(int argc,char *argv[]) {
  try { AppMain(argc,argv); }
  catch(const char *p) { cerr << p << endl; return -1; }
  catch(const string &msg) { cerr << msg.c_str() << endl; return -1; }
  catch(const exception &e) {
    cerr << "STL exception caught in main:\n" << e.what() << endl;
    return -1; }
  catch(...) { 
    cerr << "Unknown exception caught in main" << endl; return -1; }
  return 0;
}



void AppMain(int argc,char *argv[])
{
  // Process command line
  BOOM::CommandLine cmd(argc,argv,"g:CDS");
  if(cmd.numArgs()!=5)
    throw BOOM::String(
"\ncia <genezilla.iso> <alt.fasta> <projected.labels> <projected.gff> <variant-signals.gff>\n\
    <alt.fasta> = FASTA file of variant haplotype\n\
    <projected.labels> = label file produced by project-annotation\n\
    <projected.gff> = GFF file produced by project-annotation\n\
    <variant-signals.gff> = GFF file produced by find-variant-signals\n\
    <genezilla.iso> = GeneZilla configuration file, with lambda and label-matrix added\n\
  options:\n\
      -g file : dump ORF graph into file\n\
      -S : omit signal scores\n\
      -C : omit content scores\n\
      -D : omit duration scores\n\
");
  const String isochoreFilename=cmd.arg(0);
  const String fastaFilename=cmd.arg(1);
  const String labelFile=cmd.arg(2);
  const String projectedGFF=cmd.arg(3);
  const String varSigFile=cmd.arg(4);

  // Some global initialization
  alphabet=DnaAlphabet::global();
  String psaName;
  if(cmd.option('I')) psaName=cmd.optParm('I');
  bool haveEvidence=cmd.option('P');

  // Load some files
  VariantEvents variantEvents(varSigFile);
  FastaReader fastaReader(fastaFilename);
  String defline, seqString;
  fastaReader.nextSequence(defline,seqString);
  Sequence seq(seqString,alphabet);

  // Process the sequence
  float gcContent;
  int transcriptId=-1;
  EdgeFactory *edgeFactory=NULL;
  String evidenceDir;
  const int minSupport=1; // ###
  if(haveEvidence) {
    evidenceDir=cmd.optParm('P');
    edgeFactory=new FilteredEdgeFactory(NULL);
  }
  else edgeFactory=new EdgeFactory;
  CIA cia(PROGRAM_NAME,VERSION,*edgeFactory,transcriptId,variantEvents,
	  projectedGFF,labelFile);
  cia.useOneTerminusOnly();
  if(cmd.option('i')) cia.loadIsochoreBoundaries(cmd.optParm('i'));
  if(cmd.option('c')) cia.loadCpGislands(cmd.optParm('c'));
  if(cmd.option('S')) cia.omitSignalScores();
  if(cmd.option('C')) cia.omitContentScores();
  if(cmd.option('D')) cia.omitDurationScores();
  if(cmd.option('O')) cia.allowPTCs();
  EvidenceFilter *evidence=NULL;

  // Parse the substrate ID out of the defline
  BOOM::String contigId, remainder;
  BOOM::FastaReader::parseDefline(defline,contigId,remainder);

  // Predict genes and get path
  ofstream osGraph;
  bool dumpGraph=false;
  if(cmd.option('g')) {
    osGraph.open(cmd.optParm('g').c_str());
    dumpGraph=true; }
  BOOM::Stack<SignalPtr> *path=
    cia.processChunk(seq,seqString,isochoreFilename,contigId,
		     osGraph,dumpGraph,psaName);
 
  // Don't need the path; just delete it (this will also delete the
  // signal objects in it, since we use "smart pointers")
  delete path;
}



