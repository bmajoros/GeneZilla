/****************************************************************
 SpliceFeatureExtractor.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <math.h>
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
  const char *p=motif.c_str(), *q=substrate.c_str()+pos;
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
    Vector<String> fields; line.getFields(fields);
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
					       float distanceParm,
					       int maxSREdistance,
					       bool exonMotifsOnly)
  : sensor(sensor), motifs(motifs), distanceParm(distanceParm),
    consensusOffset(sensor.getConsensusOffset()),
    contextWindowLen(sensor.getContextWindowLength()),
    maxSREdistance(maxSREdistance), exonMotifsOnly(exonMotifsOnly)
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
  const int windowEnd=windowPos+contextWindowLen;
  if(windowPos<0 || windowEnd>L) return false;
  intrinsicSiteScore=sensor.getLogP(seq,seqStr,windowPos);
  regulatoryScore=computeRegScore(seqStr,consensusPos);
  return true;
}



float SpliceFeatureExtractor::computeRegScore(const String &seq,int consensusPos)
{
  const String consensus=seq.substring(consensusPos,2);
  const bool includeLeft=!exonMotifsOnly || consensus=="GT";
  const bool includeRight=!exonMotifsOnly || consensus=="AG";
  const int numMotifs=motifs.numMotifs();
  const int L=seq.length();
  float score=0;
  int begin=consensusPos-maxSREdistance, end=consensusPos+maxSREdistance;
  if(begin<0) begin=0;
  if(end>L) end=L;
  int sampleSize=0;
  for(int pos=begin ; pos<end ; ++pos) {
    if(pos<consensusPos && !includeLeft) continue;
    if(pos>consensusPos && !includeRight) continue;
    ++sampleSize;
    for(int i=0 ; i<numMotifs ; ++i) {
      const RegulatoryMotif &motif=motifs.getIthMotif(i);
      if(motif.hit(seq,pos)) {
	const int distance=abs(int(consensusPos-pos+motif.motif.length()/2));
	const float distanceScore=pow(distanceParm,distance);
	score+=motif.score*distanceScore;
	//cout<<"\t\t"<<pos<<"\t"<<i<<"\t"<<motif.motif<<"\t"<<motif.score<<"\t"<<distanceScore<<endl;
      }
    }
  }
  //return score;
  return score/sampleSize;
}


