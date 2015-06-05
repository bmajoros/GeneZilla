/****************************************************************
 GeneZilla-CIA
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_CIA_H
#include "cia.H"
#include "GZilla.H"
#include "SignalLabelingProfile.H"
#include "VariantEvents.H"
#include "ReferenceAnnotation.H"
#include "SignalStream.H"
#include "ConstraintIntervals.H"

/****************************************************************
                            class CIA
 ****************************************************************/

class CIA : public GeneZilla {
public:
  CIA(const BOOM::String &PROGRAM_NAME,const BOOM::String &VERSION,
      EdgeFactory &,int &transcriptId,VariantEvents &,const String &
      projectedGFFfile,const String &labelFile);
  virtual ~CIA();
  //virtual int main(int argc,char *argv[]);
  virtual BOOM::Stack<SignalPtr> * processChunk(const Sequence &seq,
				      const BOOM::String &seqStr,
				      const BOOM::String &isoFilename,
				      const BOOM::String &substrateId,
				      ostream &osGraph,
				      bool dumpGraph,
					String psaFilename);

protected:
  ReferenceAnnotation *refAnno;
  VariantEvents &events;
  SignalStream signalStream;
  ConstraintIntervals *constraints;
  String projectedGFF, labelFile; // filenames
  float priorWeight;
  Array1D<SignalLabelingProfile> signalLabelingProfiles;

  virtual BOOM::Stack<SignalPtr> * mainAlgorithm(const Sequence &,
						 const BOOM::String &,
						 ostream &osGraph,
						 bool dumpGraph,
						 String psaFilename);
  virtual void updateAccumulators(const Sequence &,const BOOM::String &,
				  int pos,Symbol,char);
  void scoreSignalPrior(SignalPtr);
  void initSignalLabelingProfiles();
  void buildParseGraph(const Sequence &seq,const BOOM::String &str);
};


#endif
