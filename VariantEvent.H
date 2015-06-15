/****************************************************************
 VariantEvent.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_VariantEvent_H
#define INCL_VariantEvent_H
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/Comparator.H"
#include "BOOM/Interval.H"
#include "SignalType.H"
using namespace std;
using namespace BOOM;


enum VariantSignalType {
  VAR_SIG_ATG,
  VAR_SIG_TAG,
  VAR_SIG_GT,
  VAR_SIG_AG,
  VAR_SIG_INDEL
};
VariantSignalType varSigTypeFromString(const String &);
SignalType toSignalType(VariantSignalType);
bool isSignal(VariantSignalType);

enum VariantEventType {
  VAR_EVENT_GAIN,
  VAR_EVENT_LOSS,
  VAR_EVENT_INDEL
};
VariantEventType varEventTypeFromString(const String &);

class VariantEvent {
public:
  VariantEvent(VariantSignalType,VariantEventType,int pos,int length);
  VariantSignalType getSignalType() const;
  VariantEventType getEventType() const;
  int getPosition() const;
  int begin() const;
  int end() const;
  int getLength() const;
  Interval asInterval() const { return Interval(pos,pos+length); }
  bool isGain() const;
  bool isLoss() const;
  bool isIndel() const;
  bool isStartCodon() const;
  bool isStopCodon() const;
  bool isDonor() const;
  bool isAcceptor() const;
protected:
  VariantSignalType varSignalType;
  VariantEventType varEventType;
  int pos, length;
};


struct VariantEventComparator : public Comparator<VariantEvent> {
  bool equal(VariantEvent &a,VariantEvent &b);
  bool greater(VariantEvent &a,VariantEvent &b);
  bool less(VariantEvent &a,VariantEvent &b);
};

#endif

