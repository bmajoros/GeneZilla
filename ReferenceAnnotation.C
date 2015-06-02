/****************************************************************
 ReferenceAnnotation.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "ReferenceAnnotation.H"
#include "BOOM/GffReader.H"
using namespace std;
using namespace BOOM;

ReferenceAnnotation::ReferenceAnnotation()
  : matrix(NULL), contentRegions(NULL)
{
  // ctor
}



ReferenceAnnotation::ReferenceAnnotation(const String &annoGFFfile,
					 const String &labelingFile,
					 const String &matrixFile)
{
  loadMatrix(matrixFile);
  loadLabeling(labelingFile);
  loadGFF(annoGFFfile,labeling.length());
}



ReferenceAnnotation::~ReferenceAnnotation()
{
  delete matrix;
  delete contentRegions;
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



const ContentRegions &ReferenceAnnotation::getRegions() const
{
  return *contentRegions;
}



void ReferenceAnnotation::loadGFF(const String &filename,const int seqLen)
{
  Vector<GffTranscript*> &transcripts=*GffReader::loadTranscripts(filename);
  if(transcripts.size()!=1) 
    throw "ReferenceAnnotation::loadGFF() : file must contain exactly one transcript";
  GffTranscript &transcript=*transcripts[0];
  contentRegions=new ContentRegions(transcript,seqLen);
  delete &transcripts;
}


