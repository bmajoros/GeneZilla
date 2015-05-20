/****************************************************************
 check-projection.C
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
using namespace std;
using namespace BOOM;


class Application {
public:
  Application();
  int main(int argc,char *argv[]);
private:
  GffTranscript *loadGff(const String &filename);
  void checkSpliceSites(GffTranscript &,const String &substrate);
  void checkDonor(GffExon &,const String &substrate);
  void checkAcceptor(GffExon &,const String &substrate);
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
    throw String("check-projection <ref.fasta> <ref.gff> <alt.fasta> <projected.gff>");
  const String refFasta=cmd.arg(0);
  const String refGff=cmd.arg(1);
  const String altFasta=cmd.arg(2);
  const String altGff=cmd.arg(3);

  // Load fastas
  String def, refSubstrate, altSubstrate;
  FastaReader::load(refFasta,def,refSubstrate);
  FastaReader::load(altFasta,def,altSubstrate);

  // Load GFF
  GffTranscript *refTrans=loadGff(refGff), *altTrans=loadGff(altGff);

  // Translate to proteins
  refTrans->loadSequence(refSubstrate); altTrans->loadSequence(altSubstrate);
  const String refDNA=refTrans->getSequence();
  const String altDNA=altTrans->getSequence();
  const String refProtein=ProteinTrans::translate(refDNA);
  const String altProtein=ProteinTrans::translate(altDNA);
  if(refProtein==altProtein) cout<<"proteins differ"<<endl;

  // Check for stop codons
  if(altProtein.contains("*")) cout<<"premature stop detected"<<endl;

  // Check length is divisible by 3
  if(altDNA.length()%3) cout<<"non-integral number of codons"<<endl;
  
  // Check for start codon
  if(altProtein.length()<1 || altProtein[0]!='M') cout<<"No start codon"<<endl;

  // Check splice sites
  checkSpliceSites(*altTrans,altSubstrate);

  return 0;
}



void Application::checkSpliceSites(GffTranscript &transcript,
				   const String &substrate)
{
  int numExons=transcript.getNumExons();
  for(int i=0 ; i<numExons ; ++i) {
    GffExon &exon=transcript.getIthExon(i);
    if(exon.hasDonor()) checkDonor(exon,substrate);
    if(exon.hasAcceptor()) checkAcceptor(exon,substrate);
  }
}



void Application::checkDonor(GffExon &exon,const String &substrate)
{
  if(exon.getStrand()=='+') {
    const int end=exon.getEnd();
    if(end>substrate.length()-2) return;
    String GT=substrate.substring(end,2);
    if(GT!="GT") cout<<"noncanonical donor site "<<GT<<" at "<<end<<endl;
  }
  else {
    const int begin=exon.getBegin();
    if(begin<2) return;
    String GT=ProteinTrans::reverseComplement(substrate.substring(begin-2,2));
    if(GT!="GT") cout<<"noncanonical donor site "<<GT<<" at "<<begin<<endl;
  }
}



void Application::checkAcceptor(GffExon &exon,const String &substrate)
{
  if(exon.getStrand()=='+') {
    const int begin=exon.getBegin();
    if(begin<2) return;
    String AG=substrate.substring(begin-2,2);
    if(AG!="AG") cout<<"noncanonical acceptor site "<<AG<<" at "<<begin<<endl;
  }
  else {
    const int end=exon.getEnd();
    if(end>substrate.length()-2) return;
    String AG=ProteinTrans::reverseComplement(substrate.substring(end,2));
    if(AG!="AG") cout<<"noncanonical acceptor site "<<AG<<" at "<<end<<endl;
  }
}



GffTranscript *Application::loadGff(const String &filename)
{
  GffReader reader(filename);
  Vector<GffTranscript*> *transcripts=reader.loadTranscripts();
  const int n=transcripts->size();
  if(n<1) throw filename+" contains no transcripts";
  GffTranscript *transcript=(*transcripts)[0];
  for(int i=1 ; i<n ; ++i) delete (*transcripts)[i];
  delete transcripts;
  return transcript;
}




