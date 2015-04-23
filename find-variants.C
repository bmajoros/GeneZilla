/****************************************************************
 find-variants.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/CigarString.H"
using namespace std;
using namespace BOOM;


class Application {
  int findMatch(int start,const String &ref,const String &alt,CigarString &,
		int cigarLen,int &refPos,int &altPos);
  int findMismatch(int start,const String &ref,const String &alt,CigarString &,
		   int cigarLen,int &refPos,int &altPos);
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
  if(cmd.numArgs()!=3)
    throw String("find-variants <ref.fasta> <alt.fasta> <ref-to-alt.cigar>");
  const String refFile=cmd.arg(0);
  const String altFile=cmd.arg(1);
  const String cigarFile=cmd.arg(2);

  // Load files
  CigarString rawCigar; rawCigar.load(cigarFile);
  CigarString &cigar=*rawCigar.unrollMatches();
  String refDef, altDef, refSeq, altSeq;
  FastaReader::load(refFile,refDef,refSeq);
  FastaReader::load(altFile,altDef,altSeq);
  int refLen=refSeq.length(), altLen=altSeq.length();

  // Traverse alignment
  const int L=cigar.length(); int refBegin=0, altBegin=0;
  for(int cigarPos=0 ; cigarPos<L ; ) {
    cigarPos=findMismatch(cigarPos,refSeq,altSeq,cigar,L,refBegin,altBegin);
    int refEnd=refBegin, altEnd=altBegin;
    cigarPos=findMatch(cigarPos,refSeq,altSeq,cigar,L,refEnd,altEnd);
    if(refBegin<refLen)
      cout<<"ref:\t"<<refBegin<<"\t"<<refEnd<<"\talt:\t"<<altBegin<<"\t"<<altEnd<<endl;
    refBegin=refEnd; altBegin=altEnd;
  }
  delete &cigar;

  return 0;
}



int Application::findMismatch(int cigarPos,const String &ref,const String &alt,
			   CigarString &cigar,int L,int &refPos,int &altPos)
{
  while(cigarPos<L && 
	cigar[cigarPos].type==CIGAR_MATCH && 
	ref[refPos]==alt[altPos])
    ++cigarPos, ++refPos, ++altPos;
  return cigarPos;
}



int Application::findMatch(int cigarPos,const String &ref,const String &alt,
			   CigarString &cigar,int L,int &refPos,int &altPos)
{
  while(cigarPos<L && 
	(cigar[cigarPos].type!=CIGAR_MATCH ||
	 ref[refPos]!=alt[altPos])) {
    CigarOp op=cigar[cigarPos++];
    switch(op.type) {
    case CIGAR_MATCH: ++refPos; ++altPos; break;
    case CIGAR_INSERT: altPos+=op.rep; break;
    case CIGAR_DELETE: refPos+=op.rep; break;
    }
  }
  return cigarPos;
}



