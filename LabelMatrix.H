/****************************************************************
 LabelMatrix.H
 Copyright (C)2014 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_LabelMatrix_H
#define INCL_LabelMatrix_H
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/Array2D.H"
#include "Labeling.H"
using namespace std;
using namespace BOOM;


class LabelMatrix {
public:
  LabelMatrix(const String &filename);
  float operator()(GeneModelLabel from,GeneModelLabel to);
protected:
  Array2D<float> M;
  void load(const String &filename);
};

#endif

