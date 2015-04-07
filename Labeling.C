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
  os<<labelToString(lab);
}



String labelToString(GeneModelLabel lab)
{
  switch(lab)
    {
    case LABEL_NONE:       return "?"; 
    case LABEL_INTERGENIC: return "N"; 
    case LABEL_INTRON:     return "I"; 
    case LABEL_EXON_0:     return "E0";
    case LABEL_EXON_1:     return "E1";
    case LABEL_EXON_2:     return "E2";
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



GeneModelLabel getExonLabel(int phase)
{
  switch(phase)
    {
    case 0: return LABEL_EXON_0;
    case 1: return LABEL_EXON_1;
    case 2: return LABEL_EXON_2;
    }
}



int getExonPhase(GeneModelLabel lab)
{
  switch(lab)
    {
    case LABEL_EXON_0: return 0;
    case LABEL_EXON_1: return 1;
    case LABEL_EXON_2: return 2;
    default: throw labelToString(lab)+" is not an exon";
    }
}



bool isExon(GeneModelLabel lab)
{
  switch(lab)
    {
    case LABEL_EXON_0:
    case LABEL_EXON_1:
    case LABEL_EXON_2:
      return true;
    }
  return false;
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


