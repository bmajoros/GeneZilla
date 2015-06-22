/****************************************************************
 NMD.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_NMD_H
#define INCL_NMD_H
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/GffTranscript.H"
using namespace std;
using namespace BOOM;

enum PTC_TYPE {
  PTC_NONE,
  PTC_NMD,
  PTC_TRUNCATION
};

class NMD {
public:
  NMD(int dist=50);
  PTC_TYPE predict(GffTranscript &,const String &);
protected:
  int distParm;
};

#endif

