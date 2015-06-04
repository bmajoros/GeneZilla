/****************************************************************
 ConstraintInterval.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "ConstraintInterval.H"
using namespace std;
using namespace BOOM;

ConstraintInterval::ConstraintInterval()
{
  // ctor
}



ConstraintInterval::ConstraintInterval(const Interval &interval,bool constrained)
  : interval(interval), constrained(constrained)
{
  // ctor
}



bool ConstraintInterval::isConstrained() const
{
  return constrained;
}



const Interval &ConstraintInterval::getInterval() const
{
  return interval;
}



Interval &ConstraintInterval::getInterval()
{
  return interval;
}



bool ConstraintInterval::contains(int pos) const
{
  return interval.contains(pos);
}



