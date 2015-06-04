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
}



void ConstraintIntervals::insert(ConstraintInterval I)
{
  
}



void ConstraintIntervals::reset()
{
  currentElem=0;
}



bool ConstraintIntervals::isConstrained(int pos)
{
  ConstraintInterval &interval=intervals[currentElem];
  if(interval.contains(pos)) return interval.isConstrained();
  return intervals[++currentElem].isConstrained();
}



