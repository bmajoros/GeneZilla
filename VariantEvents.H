/****************************************************************
 VariantEvents.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_VariantEvents_H
#define INCL_VariantEvents_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/Set.H"
#include "VariantEvent.H"
using namespace std;
using namespace BOOM;


class VariantEvents {
public:
  VariantEvents(const String &gffFile);
  int numEvents() const;
  VariantEvent &getIthEvent(int i);
  void eventsInInterval(const Interval &,Set<const VariantEvent*> &) const;
protected:
  Vector<VariantEvent> events;
  void load(const String &filename);
  void sort();
};

#endif

