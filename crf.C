/****************************************************************
 crf.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "CRF.H"
#include "crf.H"
#include "EdgeFactory.H"
#include "BOOM/FastaReader.H"
#include "BOOM/Constants.H"
using namespace std;

#ifdef EXPLICIT_GRAPHS
//#error Please edit the file genezilla.H, comment out the definition of EXPLICIT_GRAPHS, issue a "make clean", and recompile this project
#endif

#ifdef FORCE_SPECIFIC_SIGNALS
#error Please edit the file genezilla.H, comment out the definition of FORCE_SPECIFIC_SIGNALS, issue a "make clean", and recompile this project
#endif

static const char *PROGRAM_NAME="genezilla-CRF";
static const char *VERSION="0.1";
Alphabet alphabet;
int frame; // ### CAUTION: this is required by older code; to be removed

void AppMain(int argc,char *argv[]);

int main(int argc,char *argv[])
  {
    try
      {
	AppMain(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
	return -1;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
	return -1;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
	return -1;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
	return -1;
      }
    return 0;
  }



void AppMain(int argc,char *argv[])
{
  // Process command line
  BOOM::CommandLine cmd(argc,argv,"SCD");
  if(cmd.numArgs()!=5)
    throw BOOM::String(
    "\ncrf <*.iso> <label.matrix> <lambda> <priorlabels.txt> <*.fasta>\n\
       options:\n\
          -S : omit signal scores\n\
          -C : omit content scores\n\
          -D : omit duration scores\n\
");
  BOOM::String isochoreFilename=cmd.arg(0);
  String matrixFile=cmd.arg(1);
  float lambda=cmd.arg(2).asFloat();
  String labelFile=cmd.arg(3);
  BOOM::String fastaFilename=cmd.arg(4);

  Labeling priorLabels(labelFile);
  LabelMatrix labelMatrix(matrixFile);
  
  alphabet=DnaAlphabet::global();
  int minSeqLen=
    (cmd.option('s') ? cmd.optParm('s').asInt() : 1);
  String psaName;
  if(cmd.option('I')) psaName=cmd.optParm('I');
  bool haveEvidence=cmd.option('P');

  // Process all the contigs
  BOOM::FastaReader fastaReader(fastaFilename);
  BOOM::String defline, seqString;
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
  CRF crf(PROGRAM_NAME,VERSION,*edgeFactory,transcriptId,labelMatrix);
  if(cmd.option('i')) crf.loadIsochoreBoundaries(cmd.optParm('i'));
  if(cmd.option('c')) crf.loadCpGislands(cmd.optParm('c'));
  if(cmd.option('S')) crf.omitSignalScores();
  if(cmd.option('C')) crf.omitContentScores();
  if(cmd.option('D')) crf.omitDurationScores();
  if(cmd.option('O')) crf.allowPTCs();
  EvidenceFilter *evidence=NULL;
  while(fastaReader.nextSequence(defline,seqString))
    {
      Sequence seq(seqString,DnaAlphabet::global());
      if(seq.getLength()<minSeqLen) continue;

      // Parse the substrate ID out of the defline
      BOOM::String contigId, remainder;
      BOOM::FastaReader::parseDefline(defline,contigId,remainder);
      cerr<<"processing substrate "<<contigId<<"..."<<endl;

      // Predict genes and get path
      ofstream osGraph;
      BOOM::Stack<SignalPtr> *path=
	crf.processChunk(seq,seqString,isochoreFilename,contigId,
			       osGraph,false,psaName,priorLabels);

      // Don't need the path; just delete it (this will also delete the
      // signal objects in it, since we use "smart pointers")
      delete path;
    }
}



