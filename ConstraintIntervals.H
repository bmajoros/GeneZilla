/****************************************************************
 ConstraintIntervals.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ConstraintIntervals_H
#define INCL_ConstraintIntervals_H
#include <iostream>
#include "BOOM/Vector.H"
#include "ConstraintInterval.H"
using namespace std;
using namespace BOOM;


class ConstraintIntervals {
public:
  ConstraintIntervals(int seqLen); // creates an initial constrained interval
  void insert(ConstraintInterval); // splits an existing interval
  void reset();
  bool isConstrained(int pos);
protected:
  Vector<ConstraintInterval> intervals; // always sorted
  int currentElem;
};

#endif

