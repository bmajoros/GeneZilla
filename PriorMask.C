/****************************************************************
 PriorMask.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "PriorMask.H"
using namespace std;
using namespace BOOM;



PriorMask::PriorMask(const Interval &featureInterval)
  : interval(interval), mask(featureInterval.length())
{
  mask.setAllTo(false);
}



bool PriorMask::lookup(int geneRelativeCoord) const
{
  return mask[mapToLocal(geneRelativeCoord)];
}



void PriorMask::mask(int geneRelativeCoord)
{
  mask[mapToLocal(geneRelativeCoord)]=true;
}



int PriorMask::mapToLocal(int geneRelativeCoord) const
{
  return geneRelativeCoord-interval.getBegin();
}


