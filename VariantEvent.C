/****************************************************************
 VariantEvent.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "VariantEvent.H"
using namespace std;
using namespace BOOM;


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


