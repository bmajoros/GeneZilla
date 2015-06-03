/*
 Signal.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#include <iostream>
#include "BOOM/Constants.H"
#include "Signal.H"
#include "SignalSensor.H"
#include "SignalTypeProperties.H"
#include "Propagator.H"
#include "GarbageCollector.H"
#include "Edge.H"


Signal::Signal(int contextWindowPosition,double signalScore,
	       SignalSensor &sensor,GarbageCollector &garbageCollector)
  : contextWindowPosition(contextWindowPosition),
    sensor(sensor),
    signalScore(signalScore)
{
  predecessors[0]=predecessors[1]=predecessors[2]=NULL;

  initializePropagators(signalScore);

#ifdef EXPLICIT_GRAPHS
  garbageCollector.addSignal(this);
#endif
}



Signal::~Signal()
{
  delete propagators[0];
  delete propagators[1];

#ifdef EXPLICIT_GRAPHS
  BOOM::Set<Edge*>::iterator cur=edgesIn.begin(), end=edgesIn.end();
  for(; cur!=end ; ++cur)
    {
      Edge *edge=*cur;
      SignalPtr pred=edge->getLeft();
      pred->dropEdgeOut(edge);
      delete edge;
    }
  cur=edgesOut.begin(), end=edgesOut.end();
  for(; cur!=end ; ++cur)
    {
      Edge *edge=*cur;
      SignalPtr pred=edge->getRight();
      pred->dropEdgeIn(edge);
      delete edge;
    }
#endif
}



void Signal::dropSignalScores()
{
  signalScore=0;
}


void Signal::dropContentScores() 
{
#ifdef EXPLICIT_GRAPHS
  BOOM::Set<Edge*>::iterator cur=edgesIn.begin(), end=edgesIn.end();
  for(; cur!=end ; ++cur) {
    Edge *edge=*cur;
    edge->dropScores();
  }
#endif
}


#ifdef EXPLICIT_GRAPHS
void Signal::dropEdgeIn(Edge *edge)
{
  edgesIn.erase(edge);
}
#endif



#ifdef EXPLICIT_GRAPHS
void Signal::dropEdgeOut(Edge *edge)
{
  edgesOut.erase(edge);
}
#endif



Propagator &Signal::getPropagator(ContentType t)
{
  int index=
    SignalTypeProperties::global.whichPropagator(sensor.getSignalType(),t);
  return *propagators[index];
}



void Signal::setPredecessor(int phase,SignalPtr s)
{
  predecessors[phase]=s;
}



SignalPtr Signal::getPredecessor(int phase)
{
  return predecessors[phase];
}



int Signal::getContextWindowPosition()
{
  return contextWindowPosition;
}



int Signal::getConsensusPosition() const
{
  return contextWindowPosition+sensor.getConsensusOffset();
}



int Signal::getContextWindowLength()
{
  return sensor.getContextWindowLength();
}



int Signal::getConsensusLength()
{
  return sensor.getConsensusLength();
}



int Signal::frameOfBaseFollowingConsensus()
{
  return
    (contextWindowPosition+
     sensor.getConsensusOffset()+
     sensor.getConsensusLength())
    %3;
}



int Signal::posOfBaseFollowingConsensus()
{
  return
    contextWindowPosition+
    sensor.getConsensusOffset()+
    sensor.getConsensusLength();
}



BOOM::Set<ContentType> &Signal::belongsInWhichQueues()
{
  SignalType t=sensor.getSignalType();
  return SignalTypeProperties::global.belongsInWhichQueues(t);
}



BOOM::Set<ContentType> &Signal::linksBackToWhichQueues()
{
  SignalType t=sensor.getSignalType();
  return SignalTypeProperties::global.linksBackToWhichQueues(t);
}



Strand Signal::getStrand()
{
  return sensor.getStrand();
}



void Signal::initializePropagators(double signalScore)
{
  // Propagators start at the last base of the context window
  int propagatorPosition=
    contextWindowPosition+
    sensor.getContextWindowLength()-1;

  // Initialize the propagators to the signalScore (inductive score
  // from predecessor will have to be added later when we know what
  // that predecessor is)
  SignalType t=sensor.getSignalType();
  BOOM::Set<ContentType> &queueTypes=
    SignalTypeProperties::global.belongsInWhichQueues(t);
  int n=queueTypes.size();
  BOOM::Set<int> &phases=
    SignalTypeProperties::global.getAllowablePhases(sensor.getSignalType());
  for(int i=0 ; i<n ; ++i)
    {
      Propagator &p=*new Propagator(propagatorPosition);
      propagators[i]=&p;
      for(int i=0 ; i<3 ; ++i)
	p[i]=(phases.isMember(i) ? signalScore : NEGATIVE_INFINITY);
      //if(!isFinite(p[0]) && !isFinite(p[1]) && !isFinite(p[2])) throw "Z";
    }
  if(n<2) propagators[1]=NULL;
  //cout<<"Z "<<*this<<" "<<*propagators[0];
  //if(propagators[1]) cout<<" "<<*propagators[1]<<endl;
  //else cout<<endl;
}



SignalType Signal::getSignalType() const
{
  return sensor.getSignalType();
}



double Signal::contextWindowScore()
{
  return signalScore;
}



SignalSensor &Signal::getSensor() 
{
  return sensor;
}



double &Signal::getInductiveScore(int phase) 
{
  return inductiveScores[phase];
}



double Signal::priorInductiveScore(int phase) 
{
  return inductiveScores[phase]-signalScore;
}



double Signal::posteriorInductiveScore(int phase) 
{
  return inductiveScores[phase];
}



double Signal::precedingFeatureScore(int phase)
{
  return priorInductiveScore(phase)-
    predecessors[phase]->posteriorInductiveScore(phase);
}



#ifdef EXPLICIT_GRAPHS
void Signal::addEdgeIn(Edge *edge)
{
  edgesIn.insert(edge);
}
#endif



#ifdef EXPLICIT_GRAPHS
void Signal::addEdgeOut(Edge *edge)
{
  edgesOut.insert(edge);
}
#endif



#ifdef EXPLICIT_GRAPHS
BOOM::Set<Edge*> &Signal::getEdgesIn()
{
  return edgesIn;
}
#endif



#ifdef EXPLICIT_GRAPHS
BOOM::Set<Edge*> &Signal::getEdgesOut()
{
  return edgesOut;
}
#endif




bool Signal::isLeftTerminus()
{
  //return edgesIn.size()==0;
  return getConsensusPosition()<0;
}



#ifdef EXPLICIT_GRAPHS
bool Signal::isRightTerminus()
{
  return edgesOut.size()==0;
}
#endif



void Signal::printOn(ostream &os) const
{
  //os<<getSignalType()<<"@"<<getConsensusPosition();
  os<<getSignalType()<<"@"<<getConsensusPosition()<<":"<<*propagators[propagators[1] ? 1 :0];
}



ostream &operator<<(ostream &os,const Signal &signal)
{
  signal.printOn(os);
  return os;
}



#ifdef EXPLICIT_GRAPHS
Edge *Signal::findEdgeInFrom(SignalPtr pred)
{
  BOOM::Set<Edge*>::iterator cur=edgesIn.begin(), end=edgesIn.end();
  for(; cur!=end ; ++cur)
    {
      Edge *edge=*cur;
      if(edge->getLeft()==pred) return edge;
    }
  return NULL;
}
#endif



bool SignalPosComparator::equal(Signal* &a,Signal* &b)
{
  return a->getContextWindowPosition() == b->getContextWindowPosition();
}



bool SignalPosComparator::greater(Signal* &a,Signal* &b)
{
  return a->getContextWindowPosition() > b->getContextWindowPosition();
}


bool SignalPosComparator::less(Signal* &a,Signal* &b)
{
  return a->getContextWindowPosition() < b->getContextWindowPosition();
}

