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
  : interval(featureInterval), M(featureInterval.length())
{
  M.setAllTo(false);
}



bool PriorMask::lookup(int geneRelativeCoord) const
{
  return M[mapToLocal(geneRelativeCoord)];
}



void PriorMask::mask(int geneRelativeCoord)
{
  M[mapToLocal(geneRelativeCoord)]=true;
}



int PriorMask::mapToLocal(int geneRelativeCoord) const
{
  if(!interval.contains(geneRelativeCoord)) 
    throw String("bad coordinate in PriorMask::mapToLocal(): ")+
      geneRelativeCoord+" not in interval "+interval.getBegin()+
      ","+interval.getEnd();
  return geneRelativeCoord-interval.getBegin();
}


