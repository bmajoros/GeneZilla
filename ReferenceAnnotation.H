/****************************************************************
 ReferenceAnnotation.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ReferenceAnnotation_H
#define INCL_ReferenceAnnotation_H
#include <iostream>
#include "BOOM/Vector.H"
#include "Labeling.H"
#include "LabelMatrix.H"
#include "Signal.H"
#include "ContentRegions.H"

using namespace std;
using namespace BOOM;


class ReferenceAnnotation {
public:
  ReferenceAnnotation();
  ReferenceAnnotation(const String &annoGFFfile,const String &labelingFile,
		      const String &matrixFile);
  virtual ~ReferenceAnnotation();
  const LabelMatrix &getMatrix() const;
  const Labeling &getLabeling() const;
  const ContentRegions &getRegions() const;
protected:
  LabelMatrix *matrix;
  Labeling labeling;
  Vector<Signal> signals;
  ContentRegions *contentRegions;
  void loadMatrix(const String &filename);
  void loadLabeling(const String &filename);
  void loadGFF(const String &filename,const int seqLen);
};

#endif

