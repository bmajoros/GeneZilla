/****************************************************************
 NMD.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "NMD.H"
#include "BOOM/CodonIterator.H"
using namespace std;
using namespace BOOM;



NMD::NMD(int dist)
  : distParm(dist)
{
  // ctor
}



PTC_TYPE NMD::predict(GffTranscript &transcript,const String &substrate)
{
  const int numExons=transcript.getNumExons();
  if(numExons<1) "no exons in transcript";
  const int lastExonLen=transcript.getIthExon(numExons-1).length();
  const int lastEJC=transcript.getSplicedLength()-lastExonLen;
  CodonIterator iter(transcript,substrate);
  Codon codon;
  while(iter.nextCodon(codon))
    if(codon.isStop()) {
      const int distance=lastEJC-codon.splicedCoord;
      return (distance>=distParm) ? PTC_NMD : PTC_TRUNCATION;
    }
  return PTC_NONE;
}
