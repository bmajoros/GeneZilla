/****************************************************************
 PriorMask.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_PriorMask_H
#define INCL_PriorMask_H
#include <iostream>
#include "BOOM/Array1D.H"
#include "BOOM/Interval.H"
using namespace std;
using namespace BOOM;


class PriorMask {
public:
  PriorMask(const Interval &featureInterval);
  bool lookup(int geneRelativeCoord) const;
  void mask(int geneRelativeCoord);
protected:
  Array1D<bool> M;
  Interval interval;
  int mapToLocal(int geneRelativeCoord) const;
};

#endif

