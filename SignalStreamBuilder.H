/****************************************************************
 SignalStreamBuilder.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SignalStreamBuilder_H
#define INCL_SignalStreamBuilder_H
#include <iostream>
#include "ReferenceAnnotation.H"
#include "VariantEvents.H"
#include "ConstraintIntervals.H"
#include "SignalStream.H"
#include "Isochore.H"
using namespace std;
using namespace BOOM;


class SignalStreamBuilder {
public:
  SignalStreamBuilder(const ReferenceAnnotation &, // input
		      const VariantEvents &,       // input
		      SignalStream &,              // output
		      ConstraintIntervals &,       // output
		      Isochore *);
protected:
  const ReferenceAnnotation &refAnno;
  const VariantEvents &events;
  SignalStream &stream;
  ConstraintIntervals &constraints;
  Isochore *isochore;
  void build();
};

#endif

