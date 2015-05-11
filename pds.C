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
#include "BOOM/GffTranscriptReader.H"
using namespace std;
using namespace BOOM;


class Application
{
public:
  Application();
  int main(int argc,char *argv[]);
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
    // ctor
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=4)
    throw String("pds <ref.fasta> <ref.gff> <alt.fasta> <alt.gff>");
  const String refFasta=cmd.arg(0);
  const String refGff=cmd.arg(1);
  const String altFasta=cmd.arg(2);
  const String altGff=cmd.arg(3);

  // Load sequences
  String refDef, refSubstrate, altDef, altSubstrate;
  FastaReader::load(refFasta,refDef,refSubstrate);
  FastaReader::load(altFasta,altDef,altSubstrate);

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
    const String refSeq=refTranscript.getSequence();
    for(Vector<GffTranscript*>::iterator cur=altTranscripts.begin(), end=
	  altTranscripts.end() ; cur!=end ; ++cur) {
      GffTranscript &altTranscript=**cur;
      altTranscript.loadSequence(altSubstrate);
      const String altSeq=altTranscript.getSequence();
      sum+=compute(refSeq,altSeq);
    }
  }

  // Compute average score & print it
  float mean=sum/n;
  cout<<mean<<endl;

  // Clean up
  delete &refTranscripts; delete &altTranscripts;

  return 0;
}

