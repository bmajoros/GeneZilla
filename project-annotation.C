/****************************************************************
 project-annotation.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/CigarString.H"
#include "BOOM/GffReader.H"
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
  if(cmd.numArgs()!=5)
    throw String("project-annotation <ref.gff> <ref.fasta> <alt.fasta> <CIGAR> <out.vector>");
  const String refGffFile=cmd.arg(0);
  const String refFasta=cmd.arg(1);
  const String altFasta=cmd.arg(2);
  const String CIGAR=cmd.arg(3);const outfile=cmd.arg(4);

  CigarString cigar(CIGAR);
  GffReader gffReader(refGffFile);
  TranscriptList *transcripts=gffReader.loadTranscripts();
  int numTrans=transcripts.size();
  if(numTrans!=1) throw "Number of transcripts provided must be exactly 1";
  for(TranscriptList::iterator cur=transcripts->begin(), end=transcripts->end();
      cur!=end ; ++cur) {
    GffTranscript *transcript=*cur;

    char strand=transcript->getStrand();
    int numExons=transcript->getNumExons();

    GffExon &exon=transcript->getIthExon(i);

  }

  return 0;
}

