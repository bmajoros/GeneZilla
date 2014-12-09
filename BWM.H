/****************************************************************
 BWM.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_BWM_H
#define INCL_BWM_H
#include "WMM.H"

class BWM : public WMM
{
  int numSignificantPositions;
public:
  // This constructor is for training:
  BWM(GarbageCollector &,BOOM::Vector<TrainingSequence*> &,SignalType,
      int consensusOffset,int consensusLength,float gcContent,float alpha);
  int getNumSigPos() {return numSignificantPositions;}
};


#endif
