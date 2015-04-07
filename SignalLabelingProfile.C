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



GeneModelLabel SignalLabelingProfile::getLabel(int signalPhase,int windowPos)
{
  return M[signalPhase][windowPos];
}



void SignalLabelingProfile::init(SignalSensor &ss)
{
  const int L=ss.getContextWindowLength();
  SignalType t=ss.getSignalType();
  int offset=ss.getConsensusOffset();
  //int consLen=ss.getConsensusLength();
  M.resize(3,wlen);
  //M.setAllTo(LABEL_NONE);
  int posBaseAfterCons=offset+consLen;
  switch(t)
    {
    case ATG:       initATG(offset,L); break;
    case TAG:       initTAG(offset,L); break;
    case GT:        initGT(offset,L); break;
    case AG:        initAG(offset,L); break;
    case NEG_ATG:   initNegATG(offset,L); break;
    case NEG_TAG:   initNegTAG(offset,L); break;
    case NEG_GT:    initNegGT(offset,L); break;
    case NEG_AG:    initNegAG(offset,L); break;
    default:        M.setAllTo(LABEL_NONE); break;
    }
}



void SignalLabelingProfile::initATG(int offset,int len)
{
  M.setAllTo(LABEL_INTERGENIC);
  for(int i=offset ; i<len ; ++i)
    M[0][i]=getExonLabel((i-offset)%3);
}



void SignalLabelingProfile::initTAG(int offset,int len)
{
  M.setAllTo(LABEL_INTERGENIC);
  for(int i=offset ; i<len ; ++i)
    M[0][i]=getExonLabel(posmod(i-offset)); ###
}



void SignalLabelingProfile::initGT(int offset,int len)
{
}



void SignalLabelingProfile::initAG(int offset,int len)
{
}



void SignalLabelingProfile::initNegATG(int offset,int len)
{
}



void SignalLabelingProfile::initNegTAG(int offset,int len)
{
}



void SignalLabelingProfile::initNegGT(int offset,int len)
{
}



void SignalLabelingProfile::initNegAG(int offset,int len)
{
}




