/****************************************************************
 Labeling.C
 Copyright (C)2014 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "Labeling.H"
using namespace std;
using namespace BOOM;


Labeling::Labeling(int length)
  : A(length)
{
  // ctor
}



GeneModelLabel &Labeling::operator[](int i)
{
  return A[i];
}



Array1D<GeneModelLabel> &Labeling::asArray()
{
  return A;
}



