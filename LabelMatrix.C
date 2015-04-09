/****************************************************************
 LabelMatrix.C
 Copyright (C)2014 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "LabelMatrix.H"
#include "BOOM/File.H"
using namespace std;
using namespace BOOM;


LabelMatrix::LabelMatrix(const String &filename)
{
  load(filename);
}



float LabelMatrix::operator()(GeneModelLabel from,GeneModelLabel to)
{
  return M[from][to];
}



void LabelMatrix::load(const String &filename)
{
  M.resize(6,6);
  File f(filename);
  f.getline(); // header -- ignore
  for(int i=0 ; i<6 ; ++i) {
    String line=f.getline();
    line.trimWhitespace();
    if(line.isEmpty()) continue;
    Vector<String> &fields=*line.getFields();
    if(fields.size()<7) throw filename+" : error parsing matrix";
    for(int j=0 ; j<6 ; ++j)
      M[i][j]=fields[j+1].asFloat();
    delete &fields;
  }
}


