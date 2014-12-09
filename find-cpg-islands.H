/****************************************************************
 gilligan.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <string>
#include <iostream>
#include "BOOM/CommandLine.H"
#include "BOOM/Array1D.H"
using namespace std;

const char ZERO='0';
const char ONE='1';

class Application
{
  BOOM::String sequence, defline, substrateID;
  const char *seqStr;
  int seqLen;
  int maxSeparation;  // max separation btwn two islands to be merged
  int windowSize;     // size of sliding window
  int minIslandSize;  // minimum size of a CpG island
  float minGC;        // minimum G+C% of a CpG island
  float minRatio;     // min ratio of observed/expected CpG counts
  BOOM::Array1D<char> array;

  void fillArray();
  void eraseSmallIslands();
  void mergeCloseIslands();
  void generateOutput();
public:
  Application();
  int main(int argc,char *argv[]);
};



