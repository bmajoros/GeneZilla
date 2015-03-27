/****************************************************************
 Labeling.C
 Copyright (C)2014 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/File.H"
#include "BOOM/Map.H"
#include "Labeling.H"
using namespace std;
using namespace BOOM;


ostream &operator<<(ostream &os,GeneModelLabel lab)
{
  switch(lab)
    {
    case LABEL_NONE:       os<<"?"; break;
    case LABEL_INTERGENIC: os<<"N"; break;
    case LABEL_INTRON:     os<<"I"; break;
    case LABEL_EXON_0:     os<<"E0"; break;
    case LABEL_EXON_1:     os<<"E1"; break;
    case LABEL_EXON_2:     os<<"E2"; break;
    }
}


struct GeneModelLabelMap {
  static Map<String,GeneModelLabel> m;
  GeneModelLabelMap() {
    m["N"]=LABEL_INTERGENIC;
    m["I"]=LABEL_INTRON;
    m["E0"]=LABEL_EXON_0;
    m["E1"]=LABEL_EXON_1;
    m["E2"]=LABEL_EXON_2;
  }
};
Map<String,GeneModelLabel> GeneModelLabelMap::m;



GeneModelLabel strToLabel(const String &s)
{
  if(GeneModelLabelMap::m.isDefined(s)) return GeneModelLabelMap::m[s];
  throw s+" : unknown gene model label in strToLabel()";
}



Labeling::Labeling(int length)
  : A(length)
{
  // ctor
}



Labeling::Labeling(const String &filename)
{
  load(filename);
}



void Labeling::load(const String &filename)
{
  Vector<GeneModelLabel> V;
  File f(filename);
  while(!f.eof()) {
    String line=f.getline();
    line.trimWhitespace();
    if(line.isEmpty()) continue;
    GeneModelLabel label=strToLabel(line);
    V.push_back(label);
  }
  A=V;
}



GeneModelLabel &Labeling::operator[](int i)
{
  return A[i];
}



Array1D<GeneModelLabel> &Labeling::asArray()
{
  return A;
}



void Labeling::printOn(ostream &os) const
{
  int L=A.size();
  for(int i=0 ; i<L ; ++i) os<<A[i]<<endl;
}



ostream &operator<<(ostream &os,const Labeling &lab)
{
  lab.printOn(os);
  return os;
}



