/****************************************************************
 SignalStreamBuilder.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "SignalStreamBuilder.H"
using namespace std;
using namespace BOOM;



SignalStreamBuilder::SignalStreamBuilder(const ReferenceAnnotation &refAnno,
					 const VariantEvents &events,
					 SignalStream &stream,
					 ConstraintIntervals &constraints)
  : refAnno(refAnno), events(events), stream(stream), constraints(constraints)
{
  build();
}



void SignalStreamBuilder::build()
{
  const Vector<Signal*> &refSignals=refAnno.getSignals();
  for(Vector<Signal*>::const_iterator cur=refSignals.begin(), end=refSignals.end() ;
      cur!=end ; ++cur) stream.add(*cur);



  stream.sort();
}



