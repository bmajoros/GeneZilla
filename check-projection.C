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
#include "BOOM/CodonIterator.H"
using namespace std;
using namespace BOOM;


class Application {
public:
  Application();
  int main(int argc,char *argv[]);
private:
  Set<String> nonCanonicalGTs, nonCanonicalAGs;
  GffTranscript *loadGff(const String &filename);
  void parseNoncanonicals(const String &,Set<String> &);
  bool detectNMD(GffTranscript &altTrans,const String &altSubstrate);
  void checkSpliceSites(GffTranscript &refTrans,const String &refSubstrate,
			GffTranscript &altTrans,const String &altSubstrate);
  void checkDonor(GffExon &refExon,const String &refSubstrate,
		  GffExon &altExon,const String &altSubstrate);
  void checkAcceptor(GffExon &refExon,const String &refSubstrate,
		     GffExon &altExon,const String &altSubstrate);
  String getDonor(GffExon &,const String &substrate,int &pos);
  String getAcceptor(GffExon &,const String &substrate,int &pos);
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
  CommandLine cmd(argc,argv,"d:a:");
  if(cmd.numArgs()!=4)
    throw String("\n\
check-projection <ref.fasta> <ref.gff> <alt.fasta> <projected.gff>\n\
     -d donors = allow noncanonical donors; comma-separated list (capital letters)\n\
     -a acceptors = allow noncanonical acceptors; comma-separated list (capital letters)\n\
");
  const String refFasta=cmd.arg(0);
  const String refGff=cmd.arg(1);
  const String altFasta=cmd.arg(2);
  const String altGff=cmd.arg(3);
  if(cmd.option('d')) parseNoncanonicals(cmd.optParm('d'),nonCanonicalGTs);
  if(cmd.option('a')) parseNoncanonicals(cmd.optParm('a'),nonCanonicalAGs);

  // Load input files
  String def, refSubstrate, altSubstrate;
  FastaReader::load(refFasta,def,refSubstrate);
  FastaReader::load(altFasta,def,altSubstrate);
  GffTranscript *refTrans=loadGff(refGff), *altTrans=loadGff(altGff);

  // Translate to proteins
  refTrans->loadSequence(refSubstrate); altTrans->loadSequence(altSubstrate);
  const String refDNA=refTrans->getSequence();
  const String altDNA=altTrans->getSequence();
  const String refProtein=ProteinTrans::translate(refDNA);
  const String altProtein=ProteinTrans::translate(altDNA);
  if(refProtein!=altProtein) cout<<"proteins differ"<<endl;

  // Check for stop codons
  if(altProtein.lastChar()!='*') cout<<"missing stop codon"<<endl;
  refProtein.chop(); altProtein.chop();
  const int firstStop=altProtein.findFirst('*');
  if(firstStop>=0) {
    cout<<"premature stop at AA position "<<firstStop<<" in alt protein"<<endl;
    if(!detectNMD(*altTrans,altSubstrate))
      cout<<"truncation predicted"<<endl;
  }
  

  // Check length is divisible by 3
  if(altDNA.length()%3) cout<<"non-integral number of codons"<<endl;
  
  // Check for start codon
  if(altProtein.length()<1 || altProtein[0]!='M') cout<<"No start codon"<<endl;

  // Check splice sites
  checkSpliceSites(*refTrans,refSubstrate,*altTrans,altSubstrate);

  return 0;
}



void Application::checkSpliceSites(GffTranscript &refTranscript,
				   const String &refSubstrate,
				   GffTranscript &altTranscript,
				   const String &altSubstrate)
{
  int numExons=refTranscript.getNumExons();
  if(altTranscript.getNumExons()!=numExons)
    throw "projected transcript has different number of exons";
  for(int i=0 ; i<numExons ; ++i) {
    GffExon &refExon=refTranscript.getIthExon(i);
    GffExon &altExon=altTranscript.getIthExon(i);
    if(refExon.hasDonor()) 
      checkDonor(refExon,refSubstrate,altExon,altSubstrate);
    if(refExon.hasAcceptor()) 
      checkAcceptor(refExon,refSubstrate,altExon,altSubstrate);
  }
}



void Application::checkDonor(GffExon &refExon,const String &refSubstrate,
			     GffExon &altExon,const String &altSubstrate)
{
  int pos;
  const String refDonor=getDonor(refExon,refSubstrate,pos);
  const String altDonor=getDonor(altExon,altSubstrate,pos);
  if(altDonor==refDonor) return;
  if(altDonor=="GT") return;
  for(Set<String>::const_iterator cur=nonCanonicalGTs.begin(),
	end=nonCanonicalGTs.end() ; cur!=end ; ++cur)
    if(altDonor==*cur) return;
  cout<<"broken donor site: "<<altDonor<<" at "<<pos<<" in alt sequence"<<endl;
}



void Application::checkAcceptor(GffExon &refExon,const String &refSubstrate,
				GffExon &altExon,const String &altSubstrate)
{
  int pos;
  const String refAcceptor=getAcceptor(refExon,refSubstrate,pos);
  const String altAcceptor=getAcceptor(altExon,altSubstrate,pos);
  if(altAcceptor==refAcceptor) return;
  if(altAcceptor=="AG") return;
  for(Set<String>::const_iterator cur=nonCanonicalAGs.begin(),
	end=nonCanonicalAGs.end() ; cur!=end ; ++cur)
    if(altAcceptor==*cur) return;
  cout<<"broken acceptor site: "<<altAcceptor<<" at "<<pos<<" in alt sequence"<<endl;
}



String Application::getDonor(GffExon &exon,const String &substrate,int &pos)
{
  if(exon.getStrand()=='+') {
    const int end=exon.getEnd();
    if(end>substrate.length()-2) return "";
    return substrate.substring(pos=end,2);
  }
  else {
    const int begin=exon.getBegin();
    if(begin<2) return "";
    return ProteinTrans::reverseComplement(substrate.substring(pos=begin-2,2));
  }
}



String Application::getAcceptor(GffExon &exon,const String &substrate,int &pos)
{
  if(exon.getStrand()=='+') {
    const int begin=exon.getBegin();
    if(begin<2) return "";
    return substrate.substring(pos=begin-2,2);
  }
  else {
    const int end=exon.getEnd();
    if(end>substrate.length()-2) return "";
    return ProteinTrans::reverseComplement(substrate.substring(pos=end,2));
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



void Application::parseNoncanonicals(const String &str,Set<String> &into)
{
  Vector<String> &fields=*str.getFields(",");
  for(Vector<String>::const_iterator cur=fields.begin(), end=fields.end() ;
      cur!=end ; ++cur)
    into.insert(*cur);
  delete &fields;
}



bool Application::detectNMD(GffTranscript &transcript,
				  const String &substrate)
{
  const int numExons=transcript.getNumExons();
  if(numExons<2) return false;
  const int lastExonLen=transcript.getIthExon(numExons-1).length();
  const int lastEJC=transcript.getSplicedLength()-lastExonLen;
  CodonIterator iter(transcript,substrate);
  Codon codon;
  while(iter.nextCodon(codon))
    if(codon.isStop()) {
      const int distance=lastEJC-codon.splicedCoord;
      if(distance>=50) 
	cout<<"NMD predicted: PTC found "<<distance<<"bp from last EJC"<<endl;
      return true;
    }
  return false;
}



