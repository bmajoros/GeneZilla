/****************************************************************
 pds.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/GffReader.H"
#include "BOOM/ProteinTrans.H"
#include "BOOM/BandedSmithWaterman.H"
#include "BOOM/AminoAlphabet.H"
using namespace std;
using namespace BOOM;


class Application {
public:
  Application();
  int main(int argc,char *argv[]);
private:
  Alphabet alphabet;
  SubstitutionMatrix<float> *M;
  float gapOpen,gapExtend;
  int bandWidth;
  float compute(const Sequence &refSeq,const Sequence &altSeq);
};


int main(int argc,char *argv[])
{
  try
    {
      Application app;
      return app.main(argc,argv);
    }
  catch(const char *p)
    {
      cerr << p << endl;
    }
  catch(const string &msg)
    {
      cerr << msg.c_str() << endl;
    }
  catch(const exception &e)
    {
      cerr << "STL exception caught in main:\n" << e.what() << endl;
    }
  catch(...)
    {
      cerr << "Unknown exception caught in main" << endl;
    }
  return -1;
}



Application::Application()
{
  alphabet=static_cast<Alphabet&>(AminoAlphabet::global());
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=8)
    throw String("\n\
pds <ref.fasta> <ref.gff> <alt.fasta> <alt.gff> <matrix> <+open>\n\
    <+extend> <bandwidth_or_0>\n\
\n\
    Open and gap penalties must be positive numbers!\n\
    Higher PDS scores indicate more severe disruption (lower is better).\n\
    If bandwidth=0, non-banded alignment is performed.\n\
\n\
");
  const String refFasta=cmd.arg(0);
  const String refGff=cmd.arg(1);
  const String altFasta=cmd.arg(2);
  const String altGff=cmd.arg(3);
  String matrixFile=cmd.arg(4);
  gapOpen=-fabs(cmd.arg(5).asDouble());
  gapExtend=-fabs(cmd.arg(6).asDouble());
  int bandWidth=cmd.arg(7).asInt();
  M=new SubstitutionMatrix<float>(matrixFile,alphabet);

  // Load sequences
  String refDef, refSubstrate, altDef, altSubstrate;
  FastaReader::load(refFasta,refDef,refSubstrate);
  FastaReader::load(altFasta,altDef,altSubstrate);
  if(bandWidth<1) bandWidth=max(refFasta.length(),altFasta.length());

  // Load gene structures
  Vector<GffTranscript*> &refTranscripts=*GffReader::loadTranscripts(refGff);
  Vector<GffTranscript*> &altTranscripts=*GffReader::loadTranscripts(altGff);

  // Compute all pairwise scores
  float sum=0;
  int n=0;
  for(Vector<GffTranscript*>::iterator cur=refTranscripts.begin(), end=
	refTranscripts.end() ; cur!=end ; ++cur) {
    GffTranscript &refTranscript=**cur;
    refTranscript.loadSequence(refSubstrate);
    const Sequence refSeq(ProteinTrans::translate(refTranscript.getSequence()),alphabet);
    for(Vector<GffTranscript*>::iterator cur=altTranscripts.begin(), end=
	  altTranscripts.end() ; cur!=end ; ++cur) {
      GffTranscript &altTranscript=**cur;
      altTranscript.loadSequence(altSubstrate);
      const Sequence altSeq(ProteinTrans::translate(altTranscript.getSequence()),alphabet);
      const float score=compute(refSeq,altSeq);
      sum+=score;
      ++n;
    }
  }

  // Compute average score & print it
  float mean=sum/n;
  cout<<mean<<endl;

  // Clean up
  delete &refTranscripts; delete &altTranscripts;

  return 0;
}



float Application::compute(const Sequence &refSeq,const Sequence &altSeq)
{
  // Compute raw alignment score
  BandedSmithWaterman<float> aligner(alphabet,refSeq,altSeq,*M,gapOpen,
				      gapExtend,bandWidth);
  Alignment *alignment=aligner.fullAlignment();
  double score=alignment->getScore();
  delete alignment;

  // Normalize by maximal possible score
  float maxScore=0;
  int L=refSeq.getLength();
  for(int i=0 ; i<L ; ++i) maxScore+=(*M)(refSeq[i],refSeq[i]);
  score-=maxScore;

  return -score;
}



