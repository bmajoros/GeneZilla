/****************************************************************
 ConstraintIntervals.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "ConstraintIntervals.H"
using namespace std;
using namespace BOOM;



ConstraintIntervals::ConstraintIntervals(int seqLen)
{
  intervals.push_back(ConstraintInterval(Interval(0,seqLen),true));
  reset();
}



void ConstraintIntervals::insert(ConstraintInterval I)
{
  /*
    CASE 1: It falls completely within an existing interval.
            In that case, the existing interval may be replaced,
            split into two intervals, or truncated on either side.
    CASE 2: It overlaps at least two existing intervals.  Either of
            them may be truncated or replaced.
   */

  // First, find the first and last intervals overlapped by this new one
  const Interval interval=I.getInterval();
  const int begin=interval.getBegin(), end=interval.getEnd();
  int firstIndex=-1, secondIndex=-1, N=intervals.size();
  for(int i=0 ; i<N ; ++i) {
    const ConstraintInterval &other=intervals[i];
    if(other.contains(begin)) firstIndex=i;
    if(other.contains(end)) secondIndex=i;
  }
  if(firstIndex<0) firstIndex=0;
  if(secondIndex<0) secondIndex=N-1;

  // If the two intervals are actually the same, split that interval
  if(firstIndex==secondIndex) {
    Interval &other=intervals[firstIndex].getInterval();
    if(other==interval) intervals[firstIndex]=I; // replace
    else if(other.getBegin()==begin) {
      other.setBegin(end); // truncate
      intervals.insertByIndex(I,firstIndex);
    }
    else if(other.getEnd()==end) {
      other.setEnd(begin); // truncate
      intervals.insertByIndex(I,firstIndex+1);
    }
    else { // split into two
      Vector<ConstraintInterval> two;
      two.push_back(I); two.push_back(intervals[firstIndex]);
      other.setEnd(begin); two.back().getInterval().setBegin(end);
      intervals.insertByIndex(two,firstIndex+1);
    }
  }

  // Otherwise, perform truncation (or deletion in the boundary case)
  else {
    ConstraintInterval &firstInterval=intervals[firstIndex];
    ConstraintInterval &secondInterval=intervals[secondIndex];
    firstInterval.setEnd(begin); secondInterval.setBegin(end);
    int insertionIndex=firstIndex+1;
    int deletionIndex=firstIndex+1, deletionLen=secondIndex-firstIndex-1;
    if(secondInterval.getInterval().isEmpty()) { ++deletionLen; }
    if(firstInterval.getInterval().isEmpty()) 
      { --deletionIndex; ++deletionLen; --insertionIndex; }
    intervals.cut(deletionIndex,deletionLen);
    intervals.insertByIndex(I,insertionIndex);
  }
}



void ConstraintIntervals::reset()
{
  currentElem=0;
}



bool ConstraintIntervals::isConstrained(int pos)
{
  //  cout<<"cur="<<currentElem<<endl;
  ConstraintInterval &interval=intervals[currentElem];
  //cout<<"pos="<<pos<<" interval="<<interval.getInterval()<<endl;
  if(interval.contains(pos)) return interval.isConstrained();
  return intervals[++currentElem].isConstrained();
}



