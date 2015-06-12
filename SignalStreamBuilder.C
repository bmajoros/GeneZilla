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
					 Isochore *isochore,
					 int maxIntronScan,
					 int minExonLength)
  : refAnno(refAnno), events(events), stream(stream), constraints(constraints),
    isochore(isochore), maxIntronScan(maxIntronScan),
    minExonLength(minExonLength)
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
    if(event.isGain()) gain(event);
    else if(event.isLoss()) loss(event);
  }

  // Sort the signals
  stream.sort();

  // De-duplicate any signals that were made twice
  stream.deduplicate();
}



void SignalStreamBuilder::gain(const VariantEvent &event)
{
  const SignalType t=toSignalType(event.getSignalType());
  SignalSensor *sensor=isochore->signalTypeToSensor[t];
  const int pos=event.getPosition();
  switch(t) {
  case ATG: gainATG(pos,*sensor); break;
  case TAG: gainTAG(pos,*sensor); break;
  case GT:  gainGT(pos,*sensor);  break;
  case AG:  gainAG(pos,*sensor);  break;
  default: INTERNAL_ERROR;
  }
}



void SignalStreamBuilder::loss(const VariantEvent &event)
{
  const SignalType t=toSignalType(event.getSignalType());
  SignalSensor *sensor=isochore->signalTypeToSensor[t];
  const int pos=event.getPosition();
  switch(t) {
  case ATG: lossATG(pos,*sensor); break;
  case TAG: lossTAG(pos,*sensor); break;
  case GT:  lossGT(pos,*sensor);  break;
  case AG:  lossAG(pos,*sensor);  break;
  default: INTERNAL_ERROR;
  }
}



void SignalStreamBuilder::gainATG(int pos,SignalSensor &sensor)
{
}



void SignalStreamBuilder::gainTAG(int pos,SignalSensor &sensor)
{
}



void SignalStreamBuilder::gainGT(int pos,SignalSensor &sensor)
{
}



void SignalStreamBuilder::gainAG(int pos,SignalSensor &sensor)
{
}



void SignalStreamBuilder::lossATG(int pos,SignalSensor &sensor)
{
}



void SignalStreamBuilder::lossTAG(int pos,SignalSensor &sensor)
{
}



void SignalStreamBuilder::lossGT(int pos,SignalSensor &sensor)
{
  // First, determine interval to scan for other GT sites
  const ContentRegions &regions=refAnno.getRegions();
  const ContentRegion *prev, *next;
  if(!regions.findJunction(pos,prev,next)) INTERNAL_ERROR;
  Interval exon=prev->getInterval(), intron=next->getInterval();
  const int exonBegin=exon.getBegin(), exonEnd=exon.getEnd(),
    intronEnd=intron.getEnd();
  int begin=exonBegin+minExonLength; if(begin>exonEnd) begin=exonEnd;
  int end=exonEnd+maxIntronScan; if(end>intronEnd) end=intronEnd;

  // Perform scanning
  scan(begin,end,sensor);

  // Create unconstrained region: covers entire exon (so exon skipping
  // is an option) and a fixed portion of the following intron
  SignalSensor *AGsensor=isochore->signalTypeToSensor[AG];
  const int constraintBegin=exonBegin-2-AGsensor->getConsensusOffset();
  const Interval interval(constraintBegin,end);
  ConstraintInterval constraint(interval,false);
  constraints.insert(constraint);
}



void SignalStreamBuilder::lossAG(int pos,SignalSensor &sensor)
{
  // First, determine interval to scan for other AG sites
  const ContentRegions &regions=refAnno.getRegions();
  const ContentRegion *prev, *next;
  if(!regions.findJunction(pos+2,prev,next)) INTERNAL_ERROR;
  Interval intron=prev->getInterval(), exon=next->getInterval();
  const int intronBegin=intron.getBegin(), intronEnd=intron.getEnd(),
    exonEnd=exon.getEnd();
  int begin=intronEnd-maxIntronScan; if(begin<intronBegin) begin=intronBegin;
  int end=exonEnd-minExonLength; if(end<intronEnd) end=intronEnd;

  // Perform scanning
  scan(begin,end,sensor);

  // Create unconstrained region: covers entire exon (so exon skipping
  // is an option) and a fixed portion of the preceding intron
  SignalSensor *GTsensor=isochore->signalTypeToSensor[GT];
  const int constraintEnd=
    exonEnd-GTsensor->getConsensusOffset()+GTsensor->getContextWindowLength();
  const Interval interval(begin,constraintEnd);
  ConstraintInterval constraint(interval,false);
  constraints.insert(constraint);
}



void SignalStreamBuilder::scan(int begin,int end,SignalSensor &sensor)
{
  const Sequence &seq=refAnno.getAltSeq();
  const String &seqStr=refAnno.getAltSeqStr();
  end-=sensor.getContextWindowLength();
  for(int pos=begin ; pos<end ; ++pos) {
    SignalPtr signal=sensor.detect(seq,seqStr,pos);
    if(signal) stream.add(signal);
  }
}




