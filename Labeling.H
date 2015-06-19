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
  LABEL_NONE, // represents "none" or "unknown" or "don't care"
  LABEL_INTERGENIC,
  LABEL_INTRON,
  LABEL_EXON_0, // codon position 0
  LABEL_EXON_1, // codon position 1
  LABEL_EXON_2, // codon position 2
  LABEL_EXON    // any codon position
};
static const int NumGeneModelLabels;

ostream &operator<<(ostream &,GeneModelLabel);
GeneModelLabel strToLabel(const String &);
String labelToString(GeneModelLabel);
GeneModelLabel getExonLabel(int phase);
int getExonPhase(GeneModelLabel);
bool isExon(GeneModelLabel);

class Labeling {
public:
  Labeling(int length);
  Labeling(const String &filename);
  Labeling() {}
  void load(const String &filename);
  GeneModelLabel &operator[](int);
  Array1D<GeneModelLabel> &asArray();
  void printOn(ostream &) const;
  int length() const { return A.size(); }
private:
  Array1D<GeneModelLabel> A;
};
ostream &operator<<(ostream &,const Labeling &);

#endif

