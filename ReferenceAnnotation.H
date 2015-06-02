/****************************************************************
 ReferenceAnnotation.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ReferenceAnnotation_H
#define INCL_ReferenceAnnotation_H
#include <iostream>
#include "Labeling.H"
#include "LabelMatrix.H"

using namespace std;
using namespace BOOM;


class ReferenceAnnotation {
public:
  ReferenceAnnotation();
  virtual ~ReferenceAnnotation();
  void loadMatrix(const String &filename);
  const LabelMatrix &getMatrix() const;
  void loadLabeling(const String &filename);
  const Labeling &getLabeling() const;
protected:
  LabelMatrix *matrix;
  Labeling labeling;
};

#endif

