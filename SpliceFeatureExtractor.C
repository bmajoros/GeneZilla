/****************************************************************
 SpliceFeatureExtractor.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "SpliceFeatureExtractor.H"
using namespace std;
using namespace BOOM;



/****************************************************************
                       RegulatoryMotif
 ****************************************************************/
bool RegulatoryMotif::hit(const String &substrate,int pos) const
{
  const int substrateLen=substrate.length(), myLen=motif.length();
  const int end=pos+myLen;
  if(end>substrateLen) return false;
  const char *p=motif.c_str(), *q=substrate.c_str();
  for(int i=pos ; i<end ; ++i, ++p, ++q) if(*p!=*q) return false;
  return true;
}



/****************************************************************
                       RegulatoryMotifs
 ****************************************************************/
RegulatoryMotifs::RegulatoryMotifs(const String &filename)
{
  File f(filename);
  while(!f.eof()) {
    const String line=f.getline();
    Vector<string> fields; line.getFields(fields);
    if(fields.size()<2) continue;
    motifs.push_back(RegulatoryMotif(fields[0],fields[1].asFloat()));
  }
}



int RegulatoryMotifs::numMotifs() const
{
  return motifs.size();
}



RegulatoryMotif &RegulatoryMotifs::getIthMotif(int i)
{
  return motifs[i];
}



/****************************************************************
                     SpliceFeatureExtractor
 ****************************************************************/
SpliceFeatureExtractor::SpliceFeatureExtractor(SignalSensor &sensor,
					       const RegulatoryMotifs &motifs,
					       float distanceParm)
  : sensor(sensor), motifs(motifs), distanceParm(distanceParm),
    consensusOffset(sensor.getConsensusOffset()),
    contextWindowLen(sensor.getContextWindowLength())
{
  // ctor
}



bool SpliceFeatureExtractor::extract(const Sequence &seq,
				     const String &seqStr,int consensusPos,
				     float &intrinsicSiteScore,
				     float &regulatoryScore)
{
  const int L=seqStr.length();
  const int windowPos=consensusPos-consensusOffset;
  intrinsicSiteScore=sensor.getLogP(seq,seqStr,windowPos);
  regulatoryScore=computeRegScore(seqStr,consensusPos);
}



float SpliceFeatureExtractor::computeRegScore(const String &seq,int consensusPos)
{
  const int numMotifs=motifs.numMotifs();
  const int L=seq.length();
  float score=0;
  for(int pos=0 ; pos<L ; ++pos) {
    for(int i=0 ; i<numMotifs ; ++i) {
      const RegulatoryMotif &motif=motifs.getIthMotif(i);
      if(motif.hit(seq,pos)) {
	const int distance=abs(consensusPos-pos+motif.motif.length()/2);
	const float distanceScore=pow(distanceParm,distance);
	score+=motif.score*distanceScore;
      }
    }
  }
  return score;
}


