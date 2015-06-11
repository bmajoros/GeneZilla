/****************************************************************
 SpliceFeatureExtractor.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SpliceFeatureExtractor_H
#define INCL_SpliceFeatureExtractor_H
#include <iostream>
#include "BOOM/Vector.H"
using namespace std;
using namespace BOOM;



/****************************************************************
                       RegulatoryMotif
 ****************************************************************/
struct RegulatoryMotif {
  const String motif;
  const float score;
  RegulatoryMotif(const String &m,float score) : motif(m), score(score) {}
  bool hit(const String &,int pos) const;
};



/****************************************************************
                       RegulatoryMotifs
 ****************************************************************/
class RegulatoryMotifs {
public:
  RegulatoryMotifs(const String &filename);
  int numMotifs() const;
  RegulatoryMotif &getIthMotif(int i);
protected:
  Vector<RegulatoryMotifs> motifs;
};



/****************************************************************
                     SpliceFeatureExtractor
 ****************************************************************/
class SpliceFeatureExtractor {
public:
  SpliceFeatureExtractor(SignalSensor &,const RegulatoryMotifs &,float distanceParm);
  bool extract(const Sequence &seq,const String &seqStr,int consensusPos,
	       float &intrinsieSiteScore,float &regulatoryScore);
protected:
  SignalSensor &sensor;
  const int consensusOffset, contextWindowLen;
  const RegulatoryMotifs &motifs;
  const float distanceParm;
  float computeRegScore(const String &seq,int consensusPos);
};

#endif

