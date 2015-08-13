/****************************************************************
 ProjectionChecker.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "ProjectionChecker.H"
using namespace std;
using namespace BOOM;



ProjectionChecker::ProjectionChecker(const GffTranscript &refTrans,
				     const GffTranscript &altTrans,
				     const String &refSubstrate,
				     const String &altSubstrate,
				     const Labeling &labeling,
				     const Set<String> &nonCanonicalGTs, 
				     const Set<String> &nonCanonicalAGs)
  : refTrans(refTrans), altTrans(altTrans),
    refSubstrate(refSubstrate), altSubstrate(altSubstrate),
    labeling(labeling), nonCanonicalGTs(nonCanonicalGTs),
    nonCanonicalAGs(nonCanonicalAGs)
{
  // ctor
}



bool ProjectionChecker::checkSpliceSites(bool quiet)
{
  int numExons=refTrans.getNumExons();
  if(altTrans.getNumExons()!=numExons) {
    if(!quiet) cout<<"Projected transcript has different number of exons"<<endl;
    return false;
  }
  bool ok=true;
  for(int i=0 ; i<numExons ; ++i) {
    GffExon &refExon=refTrans.getIthExon(i);
    GffExon &altExon=altTrans.getIthExon(i);
    if(refExon.hasDonor()) 
      ok=checkDonor(refExon,refSubstrate,altExon,altSubstrate,quiet) && ok;
    if(refExon.hasAcceptor()) 
      ok=checkAcceptor(refExon,refSubstrate,altExon,altSubstrate,quiet) && ok;
  }
  return ok;
}



bool ProjectionChecker::checkDonor(GffExon &refExon,const String &refSubstrate,
				   GffExon &altExon,const String &altSubstrate,
				   bool quiet)
{
  int pos;
  const String refDonor=getDonor(refExon,refSubstrate,pos);
  const String altDonor=getDonor(altExon,altSubstrate,pos);
  if(altDonor==refDonor) return true;
  if(altDonor=="GT") return true;
  for(Set<String>::const_iterator cur=nonCanonicalGTs.begin(),
	end=nonCanonicalGTs.end() ; cur!=end ; ++cur)
    if(altDonor==*cur) return true;
  if(!quiet) 
    cout<<"broken donor site: "<<altDonor<<" at "<<pos<<" in alt sequence"<<endl;
  return false;
}



bool ProjectionChecker::checkAcceptor(GffExon &refExon,const String &refSubstrate,
				      GffExon &altExon,const String &altSubstrate,
				      bool quiet)
{
  int pos;
  const String refAcceptor=getAcceptor(refExon,refSubstrate,pos);
  const String altAcceptor=getAcceptor(altExon,altSubstrate,pos);
  if(altAcceptor==refAcceptor) return true;
  if(altAcceptor=="AG") return true;
  for(Set<String>::const_iterator cur=nonCanonicalAGs.begin(),
	end=nonCanonicalAGs.end() ; cur!=end ; ++cur)
    if(altAcceptor==*cur) return true;
  if(!quiet)
    cout<<"broken acceptor site: "<<altAcceptor<<" at "<<pos<<" in alt sequence"<<endl;
  return false;
}



String ProjectionChecker::getDonor(GffExon &exon,const String &substrate,int &pos)
{
  if(exon.getStrand()=='+') {
    const int end=exon.getEnd();
    if(end>substrate.length()-2) return "";
    return substrate.substring(pos=end,2);
  }
  else {
    const int begin=exon.getBegin();
    if(begin<2 || begin>substrate.length()) return "";
    return ProteinTrans::reverseComplement(substrate.substring(pos=begin-2,2));
  }
}



String ProjectionChecker::getAcceptor(GffExon &exon,const String &substrate,int &pos)
{
  if(exon.getStrand()=='+') {
    const int begin=exon.getBegin();
    if(begin<2 || begin>substrate.length()) return "";
    return substrate.substring(pos=begin-2,2);
  }
  else {
    const int end=exon.getEnd();
    if(end+2>substrate.length()) return "";
    return ProteinTrans::reverseComplement(substrate.substring(pos=end,2));
  }
}



bool ProjectionChecker::detectNMD(GffTranscript &transcript,
				  const String &substrate,
				  bool quiet)
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
	if(!quiet)
	  cout<<"NMD predicted: PTC found "<<distance<<"bp from last EJC"<<endl;
	return true;
      }
    }
  return false;
}



bool ProjectionChecker::checkFrameshifts(const Labeling &labeling,
					 const GffTranscript &transcript,
					 const String &substrate,
					 bool quiet)
{
  if(labeling.length()!=substrate.length()) {
    if(!quiet)
      cout<<"labeling and alt substrate have different lengths"<<endl;
    return false;
  }
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
    if(!quiet)
      cout<<"Frameshift detected: "<<phaseMismatches<<"/"<<total<<" = "
	  <<percentMismatch<<"% labeled exon bases change frame"<<endl;
    return false;
  }
  return true;
}


void ProjectionChecker::translate(GffTranscript &refTrans,
				  GffTranscript &altTrans,
				  String &refProtein,
				  String &altProtein)
{
  refTrans.loadSequence(refSubstrate); 
  String refDNA=refTrans.getSequence();
  refProtein=ProteinTrans::translate(refDNA);
  if(refProtein.lastChar()!='*') {
    refTrans.extendFinalExonBy3(); altTrans.extendFinalExonBy3();
    refTrans.loadSequence(refSubstrate); 
    refDNA=refTrans.getSequence();
    refProtein=ProteinTrans::translate(refDNA);
  }
  altTrans.loadSequence(altSubstrate);
  const String altDNA=altTrans.getSequence();
  altProtein=ProteinTrans::translate(altDNA);
}



bool ProjectionChecker::hasStartCodon(const String &protein)
{
  return protein.length()>0 && protein[0]=='M';
}



bool ProjectionChecker::hasStopCodon(const String &protein)
{
  return protein.lastChar()=='*';
}



bool ProjectionChecker::hasPTC(const String &protein,int &PTCpos)
{
  PTCpos=protein.findFirst('*');
  return PTCpos>=0 && PTCpos<protein.length()-1;
}



bool ProjectionChecker::geneIsWellFormed(const GffTranscript &trans,
					 const String &substrate,
					 const Set<String> &noncanonicalGTs,
					 const Set<String> &noncanonicalAGs,
					 bool &noStart,bool &noStop,
					 bool &PTC,bool &badSpliceSite,
					 String &msg)
{
  noStart=noStop=PTC=badSpliceSite=false;
  GffTranscript transcript=trans;
  transcript.loadSequence(substrate);
  String rna=transcript.getSequence();
  String protein=ProteinTrans::translate(rna);
  if(protein.length()<1) {
    noStart=noStop=true;
    return false;
  }
  if(protein[0]!='M') {
    msg+=String("\tExpected methionine, found ")+protein[0]+": "+protein+"\n";
    noStart=true;
  }
  const int firstStop=protein.findFirst('*');
  const int proteinLen=protein.length();
  if(firstStop>=0 && firstStop<proteinLen-1) {
    msg+=String("\tStop codon found at AA position ")+firstStop+", length="
      +proteinLen+"\n";
    PTC=true;
  }
  else if(protein.lastChar()!='*') {
    transcript.extendFinalExonBy3();
    transcript.loadSequence(substrate);
    String rna=transcript.getSequence();
    String protein=ProteinTrans::translate(rna);
    if(protein.lastChar()!='*') {
      msg+=String("\tExpected * at end of protein, found ")+protein.lastChar()
	+": "+protein+"\n";
      noStop=true;
    }
  }
  const int numExons=transcript.getNumExons();
  for(int i=0 ; i<numExons ; ++i) {
    GffExon &exon=transcript.getIthExon(i);
    if(exon.hasDonor()) 
      if(!checkDonor(exon,substrate,noncanonicalGTs,msg)) badSpliceSite=true;
    if(exon.hasAcceptor()) 
      if(!checkAcceptor(exon,substrate,noncanonicalAGs,msg)) badSpliceSite=true;
  }
  //cout<<"XXX\n"<<protein<<"\nPTC="<<PTC<<"\n"<<!(noStart || noStop || PTC || badSpliceSite)<<endl;//###
  return !(noStart || noStop || PTC || badSpliceSite);
}



bool ProjectionChecker::checkDonor(GffExon &refExon,const String &refSubstrate,
				   const Set<String> &nonCanonicalGTs,
				   String &msg)
{
  int pos;
  const String refDonor=getDonor(refExon,refSubstrate,pos);
  if(refDonor=="GT") return true;
  for(Set<String>::const_iterator cur=nonCanonicalGTs.begin(),
	end=nonCanonicalGTs.end() ; cur!=end ; ++cur)
    if(refDonor==*cur) return true;
  msg+=String("\tBroken donor site: ")+refDonor+" at position "+pos+"\n";
  return false;
}



bool ProjectionChecker::checkAcceptor(GffExon &refExon,const String &refSubstrate,
				      const Set<String> &nonCanonicalAGs,
				      String &msg)
{
  int pos;
  const String refAcceptor=getAcceptor(refExon,refSubstrate,pos);
  if(refAcceptor=="AG") return true;
  for(Set<String>::const_iterator cur=nonCanonicalAGs.begin(),
	end=nonCanonicalAGs.end() ; cur!=end ; ++cur)
    if(refAcceptor==*cur) return true;
  msg+=String("\tBroken acceptor site: ")+refAcceptor+" at position "+pos+"\n";
  return false;
}


