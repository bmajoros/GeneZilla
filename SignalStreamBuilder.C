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
					 ConstraintIntervals &constraints,
					 Isochore *isochore)
  : refAnno(refAnno), events(events), stream(stream), constraints(constraints),
    isochore(isochore)
{
  build();
}



void SignalStreamBuilder::build()
{
  // First, add the required reference signals to the stream
  const Vector<Signal*> &refSignals=refAnno.getSignals();
  for(Vector<Signal*>::const_iterator cur=refSignals.begin(), end=refSignals.end() ;
      cur!=end ; ++cur) stream.add(*cur);

  // Next, add de novo signals caused by mutations
  const int numEvents=events.numEvents();
  for(int i=0 ; i<numEvents ; ++i) {
    const VariantEvent &event=events.getIthEvent(i);
    if(event.isGain()) {
      const SignalType t=toSignalType(event.getSignalType());
      SignalSensor *sensor=isochore->signalTypeToSensor[t];
      const int contextWindowPos=event.getPosition()-sensor->getConsensusOffset();
      if(contextWindowPos<0) INTERNAL_ERROR;
      if(contextWindowPos+sensor->getContextWindowLength() > 
	 refAnno.getAltSeqStr().length()) INTERNAL_ERROR;
      const double score=sensor->getLogP(refAnno.getAltSeq(),refAnno.getAltSeqStr(),
					 contextWindowPos);
      Signal *signal=new Signal(contextWindowPos,score,*sensor,sensor->getGC());
      stream.add(signal);
      cout<<"added de novo "<<signal->getSignalType()<<endl;
    }
  }

  // Finally, do targeted signal sensing in regions around new or broken signals
  for(int i=0 ; i<numEvents ; ++i) {
    const VariantEvent &event=events.getIthEvent(i);
    if(event.isGain()) {
      const SignalType t=toSignalType(event.getSignalType());
      SignalSensor *sensor=isochore->signalTypeToSensor[t];
      const int pos=event.getPosition();

    }
    else if(event.isLoss()) {
    }
  }

  // Sort the signals
  stream.sort();

  // De-duplicate any signals that were made twice
  stream.deduplicate();
}



