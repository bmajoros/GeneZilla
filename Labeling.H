/****************************************************************
 Labeling.H
 Copyright (C)2014 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Labeling_H
#define INCL_Labeling_H
#include <iostream>
#include "BOOM/Array1D.H"
using namespace std;
using namespace BOOM;

enum GeneModelLabel {
  LABEL_NONE,
  LABEL_INTERGENIC,
  LABEL_INTRON,
  LABEL_EXON_0,
  LABEL_EXON_1,
  LABEL_EXON_2
};
ostream &operator<<(ostream &,GeneModelLabel);

class Labeling {
public:
  Labeling(int length);
  GeneModelLabel &operator[](int);
  Array1D<GeneModelLabel> &asArray();
  void printOn(ostream &) const;
  int length() const { return A.size(); }
private:
  Array1D<GeneModelLabel> A;
};
ostream &operator<<(ostream &,const Labeling &);

#endif

