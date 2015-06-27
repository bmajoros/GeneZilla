/****************************************************************
 find-non-splice-sites.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <fstream>
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/GffReader.H"
#include "BOOM/ProteinTrans.H"
#include "BOOM/CodonIterator.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/BandedSmithWaterman.H"
#include "BOOM/AminoAlphabet.H"
#include "IsochoreTable.H"
#include "Labeling.H"
#include "GCcontent.H"
using namespace std;
using namespace BOOM;

Alphabet alphabet=DnaAlphabet::global();

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
private:
  SignalSensor *GTsensor, *AGsensor;
  Sequence substrate;
  String substrateStr;
  int substrateLen;
  void processDonor(int consensusPos,int maxDistance);
  void processAcceptor(int consensusPos,int maxDistance);
};


int main(int argc,char *argv[]) {
  try {
    Application app;
    return app.main(argc,argv); }
  catch(const char *p) { cerr << p << endl; }
  catch(const string &msg) { cerr << msg.c_str() << endl; }
  catch(const exception &e) 
    {cerr << "STL exception caught in main:\n" << e.what() << endl;}
  catch(...)
    {cerr << "Unknown exception caught in main" << endl;}
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
  if(cmd.numArgs()!=7)
    throw String("\n\
find-non-splice-sites <genezilla.iso> <chr.fasta> <chr.gff> <max-proximal-distance> <intron-margin> <out-proximal-sites.fasta> <out-intronic-sites.fasta>\n\
");
  const String isoFile=cmd.arg(0);
  const String chrFasta=cmd.arg(1);
  const String gffFile=cmd.arg(2);
  const int maxDistance=cmd.arg(3).asInt();
  const int intronMargin=cmd.arg(4).asInt();
  const String proximalFile=cmd.arg(5);
  const String intronicFile=cmd.arg(6);

  // Load input files
  String def;
  FastaReader::load(chrFasta,def,substrateStr);
  substrate=Sequence(substrateStr,DnaAlphabet::global());
  substrateLen=substrateStr.length();
  GarbageIgnorer gc;
  IsochoreTable isochores(gc);
  isochores.load(isoFile);

  // Process all genes
  Vector<GffGene> &genes=*GffReader::loadGenes(gffFile);
  for(Vector<GffGene>::iterator cur=genes.begin(), end=genes.end() ;
      cur!=end ; ++cur) {
    GffGene &gene=*cur;
    GffTranscript *transcript=gene.longestTranscript();

    if(transcript->getStrand()!='+') continue; // ###
    transcript->loadSequence(substrateStr); 
    int numExons=transcript->numExons();

    // Get signal sensors
    String RNA=transcript->getSequence();
    const int begin=transcript->getBegin(), end=transcript->getEnd();
    String geneSeqStr=substrateStr.substring(begin,end-begin);
    float gcContent=GCcontent::get(geneSeqStr);
    Isochore *isochore=isochores.getIsochore(gcContent);
    GTsensor=isochore->signalTypeToSensor[GT];
    AGsensor=isochore->signalTypeToSensor[AG];

    // Iterate over splice sites
    for(int i=0 ; i<numExons ; ++i) {
      if(i>0 && i+1<numExons) {
	GffExon &exon=transcript->getIthExon(i);
	if(exon.hasDonor()) processDonor(exon.getEnd(),maxDistance);
	if(exon.hasAcceptor()) processAcceptor(exon.getBegin()-2,maxDistance);
      }
    }
  }
}



void Application::processDonor(int consensusPos,int maxDistance)
{
}



void Application::processAcceptor(int consensusPos,int maxDistance)
{
}




