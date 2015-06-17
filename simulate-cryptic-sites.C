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
  String refStr; // substrate (DNA, unspliced)
  String refRNA; // transcipt (spliced)
  Sequence refSeq;
  int refLen;
  GffTranscript *refTrans;
  int numExons;
  Labeling labeling;
  int nmd, truncations, sampleSize;
  int maxSampleSize;
  bool exonSkippingOnly, donorsOnly, acceptorsOnly;
  bool detectNMD(const GffTranscript &altTrans,const String &altSubstrate);
  String getDonor(GffExon &,const String &substrate,int &pos);
  String getAcceptor(GffExon &,const String &substrate,int &pos);
  void checkFrameshifts(const Labeling &,const GffTranscript &,
			const String &substrate);  
  void computeLabeling(TranscriptList *transcripts,Labeling &);
  void processDonor(int pos,int maxDistance,int whichExon);
  void processAcceptor(int pos,int maxDistance,int whichExon);
  void evaluate(GffTranscript &);
  bool refIsPseudogene(GffTranscript &);
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
  : nmd(0), truncations(0), sampleSize(0)
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"adem:");
  if(cmd.numArgs()!=4)
    throw String("\n\
simulate-cryptic-sites <genezilla.iso> <chr.fasta> <chr.gff> <max-distance>\n\
  -e = simulate only exon skipping\n\
  -d = simulate only changes in donor sites\n\
  -a = simulate only changes in acceptor sites\n\
  -m <N> = max sample size is N\n\
");
  const String isoFile=cmd.arg(0);
  const String refFasta=cmd.arg(1);
  const String refGff=cmd.arg(2);
  const int maxDistance=cmd.arg(3).asInt();
  exonSkippingOnly=cmd.option('e');
  donorsOnly=cmd.option('d');
  acceptorsOnly=cmd.option('a');
  maxSampleSize=cmd.option('m') ? cmd.optParm('m').asInt() : 0;

  // Load input files
  String def;
  FastaReader::load(refFasta,def,refStr);
  refSeq=Sequence(refStr,DnaAlphabet::global());
  refLen=refStr.length();
  GarbageIgnorer gc;
  IsochoreTable isochores(gc);
  isochores.load(isoFile);

  // Process all genes
  Vector<GffGene> &genes=*GffReader::loadGenes(refGff);
  for(Vector<GffGene>::iterator cur=genes.begin(), end=genes.end() ;
      cur!=end ; ++cur) {
    GffGene &gene=*cur;
    refTrans=gene.longestTranscript();
    if(maxSampleSize>0 && sampleSize>maxSampleSize) break;

    if(refTrans->getStrand()!='+') continue; // ###
    if(refIsPseudogene(*refTrans)) continue;
    refTrans->loadSequence(refStr); 
    numExons=refTrans->numExons();
    labeling=Labeling(refLen);

    // Get signal sensors
    refRNA=refTrans->getSequence();
    float gcContent=GCcontent::get(refRNA);
    Isochore *isochore=isochores.getIsochore(gcContent);
    GTsensor=isochore->signalTypeToSensor[GT];
    AGsensor=isochore->signalTypeToSensor[AG];

    // Iterate over splice sites
    for(int i=0 ; i<numExons ; ++i) {
      if(exonSkippingOnly) {
	if(i>0 && i+1<numExons) {
	  GffTranscript altTrans(*refTrans);
	  altTrans.deleteIthExon(i);
	  //	  cout<<altTrans.getIthExon(i).length()%3<<endl;
	  //if(altTrans.getIthExon(i).length()%3) continue;
	  {//###
	    /*
	    GffExon &prev=altTrans.getIthExon(i-1), &next=altTrans.getIthExon(i+1);
	    cout<<prev.getSequence().substr(prev.getSequence().length()-2)<<" "
		<<next.getSequence().substr(0,2)<<endl;
	    */
	  }//###
	  evaluate(altTrans);
	}
	continue;
      }
      GffExon &exon=refTrans->getIthExon(i);
      if(exon.hasDonor()) processDonor(exon.getEnd(),maxDistance,i);
      if(exon.hasAcceptor()) processAcceptor(exon.getBegin()-2,maxDistance,i);
    }
    if(sampleSize>0) {
      const int stops=nmd+truncations;
      const float percentStops=stops/float(sampleSize);
      const float nmdOverStops=nmd/float(stops);
      const float truncOverStops=truncations/float(stops);
      const float nmdOverAll=nmd/float(sampleSize);
      cout<<"NMD="<<nmd<<" trunc="<<truncations<<" sample="<<sampleSize
	  <<" %stops="<<percentStops<<" #nmd/#stops="<<nmdOverStops
	  <<" #trunc/#stops="<<truncOverStops
	  <<" #nmd/#sample="<<nmdOverAll<<endl;
    }
  }

  // Report statistics;
    if(sampleSize>0) {
      const int stops=nmd+truncations;
      const float percentStops=stops/float(sampleSize);
      const float nmdOverStops=nmd/float(stops);
      const float truncOverStops=truncations/float(stops);
      const float nmdOverAll=nmd/float(sampleSize);
      cout<<"NMD="<<nmd<<" trunc="<<truncations<<" sample="<<sampleSize
	  <<" %stops="<<percentStops<<" #nmd/#stops="<<nmdOverStops
	  <<" #trunc/#stops="<<truncOverStops
	  <<" #nmd/#sample="<<nmdOverAll<<endl;
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
	//cout<<"NMD predicted: PTC found "<<distance<<"bp from last EJC"<<endl;
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
  if(acceptorsOnly) return;
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
	GffExon &exon=altTrans.getIthExon(whichExon);
	exon.setEnd(pos);
	if(exon.length()>0) evaluate(altTrans);
      }
    }
  }
}



void Application::processAcceptor(int refPos,int maxDistance,int whichExon)
{
  if(donorsOnly) return;
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
	GffExon &exon=altTrans.getIthExon(whichExon);
	exon.setBegin(pos+2);
	if(exon.length()>0) evaluate(altTrans);
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



bool Application::refIsPseudogene(GffTranscript &trans)
{
  trans.loadSequence(refStr); 
  const String altRNA=trans.getSequence();
  const String altProtein=ProteinTrans::translate(altRNA);
  altProtein.chop();
  const int firstStop=altProtein.findFirst('*');
  return firstStop>=0;
}



void Application::evaluate(GffTranscript &trans)
{
  trans.loadSequence(refStr); 
  const String altRNA=trans.getSequence();
  const String refProtein=ProteinTrans::translate(refRNA);
  const String altProtein=ProteinTrans::translate(altRNA);
  refProtein.chop(); altProtein.chop();
  if(refProtein==altProtein) throw "identical"; // ### debugging
  const int firstStop=altProtein.findFirst('*');
  if(firstStop>=0) {
    if(maxSampleSize<1 || sampleSize<=maxSampleSize) {
      if(detectNMD(trans,refStr)) ++nmd;
      else ++truncations;
    }
  }
  ++sampleSize;
}



