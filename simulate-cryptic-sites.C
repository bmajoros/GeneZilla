/****************************************************************
 simulate-cryptic-sites.C
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
#include "BOOM/DnaAlphabet.H"
#include "IsochoreTable.H"
#include "Labeling.H"
#include "GCcontent.H"
using namespace std;
using namespace BOOM;

/*
  Statistics from DBASS:
  ===> 95% of splicing changes in DBASS activate a site less than 70bp away
  ===> 93% are less than 50bp
  ===> 90% are less than 30bp
 */

Alphabet alphabet=DnaAlphabet::global();

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
private:
  SignalSensor *GTsensor, *AGsensor;
  String refStr;
  Sequence refSeq;
  int refLen;
  GffTranscript *refTrans;
  int numExons;
  Labeling labeling;
  int nmd, nmdSampleSize;
  bool detectNMD(const GffTranscript &altTrans,const String &altSubstrate);
  //void checkSpliceSites(GffTranscript &refTrans,const String &refSubstrate,
  //			GffTranscript &altTrans,const String &altSubstrate);
  //void checkDonor(GffExon &refExon,const String &refSubstrate,
  //		  GffExon &altExon,const String &altSubstrate);
  //void checkAcceptor(GffExon &refExon,const String &refSubstrate,
  //		     GffExon &altExon,const String &altSubstrate);
  String getDonor(GffExon &,const String &substrate,int &pos);
  String getAcceptor(GffExon &,const String &substrate,int &pos);
  void checkFrameshifts(const Labeling &,const GffTranscript &,
			const String &substrate);  
  void computeLabeling(TranscriptList *transcripts,Labeling &);
  void processDonor(int pos,int maxDistance,int whichExon);
  void processAcceptor(int pos,int maxDistance,int whichExon);
  void evaluate(const GffTranscript &);
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
  : nmd(0), nmdSampleSize(0)
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=4)
    throw String("\n\
simulate-cryptic-sites <genezilla.iso> <ref.fasta> <ref.gff> <max-distance>\n\
");
  const String isoFile=cmd.arg(0);
  const String refFasta=cmd.arg(1);
  const String refGff=cmd.arg(2);
  const int maxDistance=cmd.arg(3).asInt();

  // Load input files
  String def;
  FastaReader::load(refFasta,def,refStr);
  refSeq=Sequence(refStr,DnaAlphabet::global());
  refLen=refStr.length();
  refTrans=GffReader::longestTranscript(refGff);
  if(refTrans->getStrand()!='+') 
    throw refGff+": only the forward strand is currently supported";
  refTrans->loadSequence(refStr); 
  numExons=refTrans->numExons();
  labeling=Labeling(refLen);

  // Load signal sensors
  GarbageIgnorer gc;
  IsochoreTable isochores(gc);
  isochores.load(isoFile);
  float gcContent=GCcontent::get(refStr);
  Isochore *isochore=isochores.getIsochore(gcContent);
  GTsensor=isochore->signalTypeToSensor[GT];
  AGsensor=isochore->signalTypeToSensor[AG];

  // Iterate over splice sites
  for(int i=0 ; i<numExons ; ++i) {
    GffExon &exon=refTrans->getIthExon(i);
    if(exon.hasDonor()) processDonor(exon.getEnd(),maxDistance,i);
    if(exon.hasAcceptor()) processAcceptor(exon.getBegin()-2,maxDistance,i);
  }
}



/*
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
*/



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



bool Application::detectNMD(const GffTranscript &transcript,
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
      if(distance>=50) {
	cout<<"NMD predicted: PTC found "<<distance<<"bp from last EJC"<<endl;
	return true;
      }
    }
  return false;
}



void Application::checkFrameshifts(const Labeling &labeling,
				   const GffTranscript &transcript,
				   const String &substrate)
{
  if(labeling.length()!=substrate.length()) 
    throw "labeling and alt substrate have different lengths";
  const int numExons=transcript.numExons();
  int phase=0, phaseMatches=0, phaseMismatches=0;
  for(int i=0 ; i<numExons ; ++i) {
    const GffExon &exon=transcript.getIthExon(i);
    const int begin=exon.getBegin(), end=exon.getEnd();
    for(int pos=begin ; pos<end ; ++pos) {
      const GeneModelLabel label=labeling[pos];
      if(isExon(label))
	if(phase==getExonPhase(label)) ++phaseMatches;
	else ++phaseMismatches;
      phase=(phase+1)%3;
    }
  }
  if(phaseMismatches>0) {
    const int total=phaseMismatches+phaseMatches;
    float percentMismatch=int(1000*phaseMismatches/float(total)+5/9.0)/10.0;
    cout<<"frameshift detected: "<<phaseMismatches<<"/"<<total<<" = "
	<<percentMismatch<<"% labeled exon bases change frame"<<endl;
  }
}



void Application::computeLabeling(TranscriptList *transcripts,
				  Labeling &refLab)
{
  int numTrans=transcripts->size();
  if(numTrans!=1) throw "Number of transcripts provided must be exactly 1";
  GffTranscript *transcript=(*transcripts)[0];
  int begin=transcript->getBegin(), end=transcript->getEnd();
  char strand=transcript->getStrand();
  if(strand!='+') throw "only forward-strand features are currently supported";
  int numExons=transcript->getNumExons();
  refLab.asArray().setAllTo(LABEL_INTERGENIC);
  for(int i=begin ; i<end ; ++i) refLab[i]=LABEL_INTRON;
  int phase=0;
  for(int i=0 ; i<numExons ; ++i) {
    GffExon &exon=transcript->getIthExon(i);
    begin=exon.getBegin(); end=exon.getEnd();
    for(int j=begin ; j<end ; ++j) {
      refLab[j]=getExonLabel(phase);
      phase=(phase+1)%3;
    }
  }
}



void Application::processDonor(int refPos,int maxDistance,int whichExon)
{
  const int consensusOffset=GTsensor->getConsensusOffset();
  const int contextWindowLen=GTsensor->getContextWindowLength();
  const float threshold=GTsensor->getCutoff();
  int firstPos=refPos-maxDistance; if(firstPos<0) firstPos=0;
  int lastPos=refPos+2+maxDistance; if(lastPos>refLen-2) lastPos=refLen-2;
  for(int pos=firstPos ; pos<=lastPos ; ++pos) {
    if(pos==refPos) continue;
    if(GTsensor->consensusOccursAt(refStr,pos)) {
      const int begin=pos-consensusOffset;
      const int end=begin+contextWindowLen;
      if(end>refLen) continue;
      const float logP=GTsensor->getLogP(refSeq,refStr,pos-consensusOffset);
      if(logP>=threshold) {
	GffTranscript altTrans(*refTrans);
	altTrans.getIthExon(whichExon).setEnd(pos);
	evaluate(altTrans);
      }
    }
  }
}



void Application::processAcceptor(int refPos,int maxDistance,int whichExon)
{
  const int consensusOffset=AGsensor->getConsensusOffset();
  const int contextWindowLen=AGsensor->getContextWindowLength();
  const float threshold=AGsensor->getCutoff();
  int firstPos=refPos-maxDistance; if(firstPos<0) firstPos=0;
  int lastPos=refPos+2+maxDistance; if(lastPos>refLen-2) lastPos=refLen-2;
  for(int pos=firstPos ; pos<=lastPos ; ++pos) {
    if(pos==refPos) continue;
    if(AGsensor->consensusOccursAt(refStr,pos)) {
      const int begin=pos-consensusOffset;
      const int end=begin+contextWindowLen;
      if(end>refLen) continue;
      const float logP=AGsensor->getLogP(refSeq,refStr,pos-consensusOffset);
      if(logP>=threshold) {
	GffTranscript altTrans(*refTrans);
	altTrans.getIthExon(whichExon).setBegin(pos);
	evaluate(altTrans);
      }
    }
  }
}



/*
void Application::checkMutation()
{
  // Translate to proteins
  refTrans->loadSequence(refSubstrate); 
  String refDNA=refTrans->getSequence();
  String refProtein=ProteinTrans::translate(refDNA);
  if(refProtein.lastChar()!='*') {
    cout<<"extending"<<endl;
    refTrans->extendFinalExonBy3(); altTrans->extendFinalExonBy3();
    refTrans->loadSequence(refSubstrate); 
    refDNA=refTrans->getSequence();
    refProtein=ProteinTrans::translate(refDNA);
  }
  altTrans->loadSequence(altSubstrate);
  const String altDNA=altTrans->getSequence();
  const String altProtein=ProteinTrans::translate(altDNA);
  if(refProtein!=altProtein) cout<<"proteins differ"<<endl;

  // Check for frameshifts
  checkFrameshifts(labeling,*altTrans,altSubstrate);

  // Check for stop codons
  const bool stopPresent=altProtein.lastChar()=='*';
  refProtein.chop(); altProtein.chop();
  const int firstStop=altProtein.findFirst('*');
  if(firstStop>=0) {
    cout<<"premature stop at AA position "<<firstStop<<" in alt protein"<<endl;
    if(!detectNMD(*altTrans,altSubstrate))
      cout<<"truncation predicted"<<endl;
  }
  else if(!stopPresent) cout<<"missing stop codon"<<endl; 
  
  // Check length is divisible by 3
  if(altDNA.length()%3) cout<<"non-integral number of codons"<<endl;
  
  // Check for start codon
  if(altProtein.length()<1 || altProtein[0]!='M') cout<<"No start codon"<<endl;

  // Check splice sites
  checkSpliceSites(*refTrans,refSubstrate,*altTrans,altSubstrate);

  return 0;
}
*/



void Application::evaluate(const GffTranscript &trans)
{
  trans.loadSequence(refStr); 
  const String refRNA=refTrans->getSequence();
  const String altRNA=trans.getSequence();
  const String refProtein=ProteinTrans::translate(refRNA);
  const String altProtein=ProteinTrans::translate(refRNA);
  refProtein.chop(); altProtein.chop();
  const int firstStop=altProtein.findFirst('*');
  if(firstStop>=0) {
    if(!detectNMD(trans,refStr)) ++nmd;
    ++nmdSampleSize;
  }
}



