/****************************************************************
 VariantEvents.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "VariantEvents.H"
using namespace std;
using namespace BOOM;

VariantEvents::VariantEvents(const String &filename)
{
  load(filename);
}


int VariantEvents::numEvents() const
{
  return events.size();
}



VariantEvent &VariantEvents::getIthEvent(int i)
{
  return events[i];
}



void VariantEvents::load(const String &filename)
{
  GffReader reader(filename);
  Vector<GffFeature*> *features=reader.loadFeatures();
  for(Vector<GffFeature*>::iterator cur=features->begin(), end=features->end() ; 
      cur!=end ; ++cur) {
    const GffFeature &feature=**cur;
    const VariantSignalType varSigType=varSigTypeFromString(feature.getFeatureType());
    const VariantEventType varEventType=varEventTypeFromString(feature.getSource());
    const int pos=feature.getBegin();
    events.push_back(VariantEvent(varSigType,varEventType,pos));
    delete &feature;
  }
}

