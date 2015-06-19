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
  Set<Signal*> newSignals;

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
  void enforceConstraints(SignalPtr);
  void purgeQueues();
  void reweightGraph();
  void reweight(Edge &);
  void makePriorMask(Array1D<bool> &mask,const Edge &,bool leftIsNew,
		     bool rightIsNew,const Set<const VariantEvent*> 
		     &coveredEvents);
  void maskLeft(Array1D<bool> &,const Edge &);
  void maskLeftGT(Array1D<bool> &mask,Signal *,const Edge &);
  void maskLeftAG(Array1D<bool> &mask,Signal *,const Edge &);
  void maskLeftATG(Array1D<bool> &mask,Signal *,const Edge &);
  void maskRight(Array1D<bool> &,const Edge &);
  void maskRightGT(Array1D<bool> &mask,Signal *,const Edge &);
  void maskRightAG(Array1D<bool> &mask,Signal *,const Edge &);
  void maskEvents(Array1D<bool> &,const Edge &,
		  const Set<const VariantEvent*> &coveredEvents);
  void regionOverlapping(int pos,ContentType &,Interval &);
  void mask(const Interval &,Array1D<bool> &,int maskOffset);
  void computePrior(Edge &,const Array1D<bool> &mask);
  float computePrior(const Labeling &refLabels,const Labeling &altLabels,
		     int offset,const Array1D<bool> &mask);
};


#endif
