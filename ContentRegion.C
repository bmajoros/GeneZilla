/****************************************************************
 ContentRegion.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "ContentRegion.H"
using namespace std;
using namespace BOOM;

ContentRegion::ContentRegion()
{
  // ctor
}



ContentRegion::ContentRegion(ContentType t,int begin,int end)
  : type(t), interval(begin,end)
{
  // ctor
}



ContentType ContentRegion::getType() const
{
  return type;
}



const Interval &ContentRegion::getInterval() const
{
  return interval;
}


