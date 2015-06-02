/****************************************************************
 ReferenceAnnotation.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "ReferenceAnnotation.H"
using namespace std;
using namespace BOOM;

ReferenceAnnotation::ReferenceAnnotation()
  : matrix(NULL)
{
  // ctor
}



ReferenceAnnotation::~ReferenceAnnotation()
{
  delete matrix;
}



void ReferenceAnnotation::loadMatrix(const String &filename)
{
  matrix=new LabelMatrix(filename);
}



const LabelMatrix &ReferenceAnnotation::getMatrix() const
{
  return *matrix;
}



void ReferenceAnnotation::loadLabeling(const String &filename)
{
  labeling.load(filename);
}



const Labeling &ReferenceAnnotation::getLabeling() const
{
  return labeling;
}




