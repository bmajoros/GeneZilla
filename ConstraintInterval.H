/****************************************************************
 ConstraintInterval.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ConstraintInterval_H
#define INCL_ConstraintInterval_H
#include <iostream>
#include "BOOM/Interval.H"
using namespace std;
using namespace BOOM;


class ConstraintInterval {
public:
  ConstraintInterval();
  ConstraintInterval(const Interval &,bool constrained);
  bool isConstrained() const;
  bool contains(int pos) const;
  const Interval &getInterval() const;
protected:
  bool constrained;
  Interval interval;
};

#endif

