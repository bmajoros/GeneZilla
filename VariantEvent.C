/****************************************************************
 VariantEvent.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/Exceptions.H"
#include "VariantEvent.H"
using namespace std;


VariantEvent::VariantEvent(VariantSignalType vst,VariantEventType vet,
			   int pos,int length)
  : varSignalType(vst), varEventType(vet), pos(pos), length(length)
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



bool VariantEvent::isGain() const
{
  return varEventType==VAR_EVENT_GAIN;
}



bool VariantEvent::isLoss() const
{
  return varEventType==VAR_EVENT_LOSS;
}



bool VariantEvent::isIndel() const
{
  return varEventType==VAR_EVENT_INDEL;
}



bool VariantEvent::isStartCodon() const
{
  return varSignalType==VAR_SIG_ATG;
}



bool VariantEvent::isStopCodon() const
{
  return varSignalType==VAR_SIG_TAG;
}



bool VariantEvent::isDonor() const
{
  return varSignalType==VAR_SIG_GT;
}



bool VariantEvent::isAcceptor() const
{
  return varSignalType==VAR_SIG_AG;
}



VariantSignalType varSigTypeFromString(const String &s)
{
  if(s=="start-codon") return VAR_SIG_ATG;
  if(s=="stop-codon") return VAR_SIG_TAG;
  if(s=="donor-site" || s=="donor") return VAR_SIG_GT;
  if(s=="acceptor-site" || s=="acceptor") return VAR_SIG_AG;
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



SignalType toSignalType(VariantSignalType t)
{
  switch(t) {
  case VAR_SIG_ATG:    return ATG;
  case VAR_SIG_TAG:    return TAG;
  case VAR_SIG_GT:     return GT;
  case VAR_SIG_AG:     return AG;
  case VAR_SIG_INDEL: 
  default:
    INTERNAL_ERROR;
  }
}



bool isSignal(VariantSignalType t)
{
  return t!=VAR_SIG_INDEL;
}



int VariantEvent::begin() const
{
  return pos;
}



int VariantEvent::end() const
{
  return pos+length;
}



int VariantEvent::getLength() const
{
  return length;
}



bool VariantEventComparator::equal(VariantEvent &a,VariantEvent &b)
{
  return a.begin()==b.begin();
}



bool VariantEventComparator::greater(VariantEvent &a,VariantEvent &b)
{
  return a.begin()>b.begin();
}



bool VariantEventComparator::less(VariantEvent &a,VariantEvent &b)
{
  return a.begin()<b.begin();
}



