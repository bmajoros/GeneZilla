/****************************************************************
 SignalStream.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "SignalStream.H"
#include "BOOM/VectorSorter.H"
using namespace std;
using namespace BOOM;



SignalStream::SignalStream()
  : currentIndex(0)
{
  // ctor
}



SignalStream::~SignalStream()
{
  for(Vector<Signal*>::iterator cur=signals.begin(), end=signals.end() ;
      cur!=end ; ++cur)
    delete *cur;
}



void SignalStream::add(Signal *s)
{
  signals.push_back(s);
}



void SignalStream::sort()
{
  SignalPosComparator cmp;
  VectorSorter<Signal*> sorter(signals,cmp);
  sorter.sortAscendInPlace();
}



Signal *SignalStream::detect(int position)
{
  int N=signals.size();
  if(currentIndex>=N) return NULL;
  const sigPos=signals[currentIndex]->getContextWindowPosition();
  if(position==sigPos) {
    ++currentIndex;
    return signals[currentIndex];
  }
  if(position>sigPos) INTERNAL_ERROR; // ### DEBUGGING
  return NULL;
}



void SignalStream::reset()
{
  currentIndex=0;
}

