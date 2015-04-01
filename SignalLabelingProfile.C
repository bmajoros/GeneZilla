/****************************************************************
 SignalLabelingProfile.C
 Copyright (C)2014 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "SignalLabelingProfile.H"
using namespace std;
using namespace BOOM;

SignalLabelingProfile::SignalLabelingProfile(SignalSensor &ss)
{
  init();
}



GeneModelLabel SignalLabelingProfile:getLabel(int signalPhase,int windowPos)
{
  return M[signalPhase][windowPos];
}



void SignalLabelingProfile:init(SignalSensor &ss)
{
  const int L=ss.getContextWindowLength();
  SignalType t=ss.getSignalType();
  int offset=ss.getConsensusOffset();
  int consLen=ss.getConsensusLength();
  M.resize(3,wlen);
  M.setAllTo(LABEL_NONE);
  int posBaseAfterCons=offset+consLen;
  switch(t)
    {
    case ATG:
    case TAG:
    case GT:
    case AG:
    case PROM:
    case POLYA:
    case NEG_ATG:
    case NEG_TAG:
    case NEG_GT:
    case NEG_AG:
    case NEG_PROM:
    case NEG_POLYA:
    }

}

