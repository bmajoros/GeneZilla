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




SpliceFeatureExtractor::SpliceFeatureExtractor()
{
  // ctor
}


