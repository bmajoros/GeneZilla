/****************************************************************
 VariantEvent.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "VariantEvent.H"
using namespace std;


VariantEvent::VariantEvent(VariantSignalType vst,VariantEventType vet,int pos)
  : varSignalType(vst), varEventType(vet), pos(pos)
{
  // ctor
}



VariantSignalType VariantEvent::getSignalType() const
{
  return varSignalType;
}



VariantEventType VariantEvent::getEventType() const
{
  return varEventType;
}



int VariantEvent::getPosition() const
{
  return pos;
}



VariantSignalType varSigTypeFromString(const String &s)
{
  if(s=="start-codon") return VAR_SIG_ATG;
  if(s=="stop-codon") return VAR_SIG_TAG;
  if(s=="donor-site") return VAR_SIG_GT;
  if(s=="acceptor-site") return VAR_SIG_AG;
  if(s=="deletion" || s=="insertion") return VAR_SIG_INDEL;
  throw s+" : unknown VariantSignalType";
};



VariantEventType varEventTypeFromString(const String &s)
{
  if(s=="new") return VAR_EVENT_GAIN;
  if(s=="broken") return VAR_EVENT_LOSS;
  if(s=="indel") return VAR_EVENT_INDEL;
  throw s+" : unknown VariantEventType";
};


