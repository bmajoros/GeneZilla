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
  for(int i=offset, phase=0 ; i<len ; ++i, phase=(phase+1)%3)
    M[0][i]=getExonLabel(phase);
}



void SignalLabelingProfile::initTAG(int offset,int len)
{
  M.setAllTo(LABEL_INTERGENIC);
  for(int i=offset+2, phase=2 ; i>=0 ; --i, phase=posmod(phase-1))
    M[0][i]=getExonLabel(phase);
}



void SignalLabelingProfile::initGT(int offset,int len)
{
  M.setAllTo(LABEL_INTRON);
  for(int signalPhase=0 ; signalPhase<3 ; ++signalPhase)
    for(int i=offset-1, phase=posmod(signalPhase-1) ; i>=0 ; 
	--i, phase=posmod(phase-1))
      M[signalPhase][i]=getExonLabel(phase);
}



void SignalLabelingProfile::initAG(int offset,int len)
{
  M.setAllTo(LABEL_INTRON);
  for(int signalPhase=0 ; signalPhase<3 ; ++signalPhase)
    for(int i=offset+2, phase=signalPhase ; i<len ; ++i, phase=(phase+1)%3)
      M[signalPhase][i]=getExonLabel(phase);
}



void SignalLabelingProfile::initNegATG(int offset,int len)
{
  M.setAllTo(LABEL_INTERGENIC);
  for(int i=offset+2, phase=0 ; i>=0 ; --i, phase=(phase+1)%3)
    M[0][i]=getExonLabel(phase);
}



void SignalLabelingProfile::initNegTAG(int offset,int len)
{
  M.setAllTo(LABEL_INTERGENIC);
  for(int i=offset, phase=2 ; i<len ; ++i, phase=posmod(phase-1);
    M[0][i]=getExonLabel(phase);
}



void SignalLabelingProfile::initNegGT(int offset,int len)
{
  M.setAllTo(LABEL_INTRON);
  for(int signalPhase=0 ; signalPhase<3 ; ++signalPhase)
    for(int i=offset+2, phase=signalPhase ; i<len ; ++i, phase=posmod(phase-1);
	M[signalPhase][i]=getExonLabel(phase);
}



void SignalLabelingProfile::initNegAG(int offset,int len)
{
  M.setAllTo(LABEL_INTRON);
  for(int signalPhase=0 ; signalPhase<3 ; ++signalPhase)
    for(int i=offset-1, phase=signalPhase ; i>=0 ; --i, phase=posmod(phase-1);
	M[signalPhase][i]=getExonLabel(phase);
}


/*

----ATG----   => signal phase is 0
    0120120

----TAG----   => signal phase is 0
2012012

----GT----    => signal phase is 1
0120

----AG----    => signal phase is 2
      2012

----CAT----   => signal phase is 2 (reverse-strand start codon)
0210210

----CTA----   => signal phase is 2 (reverse-strand stop codon)
    2102102

----AC----    => signal phase is 0 (reverse-strand donor)
      0210

----CT----    => signal phase is 1 (reverse-strand acceptor)
2102

 */


