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
#include "PriorMask.H"

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
  void makePriorMask(PriorMask &,const Edge &,bool leftIsNew,
		     bool rightIsNew,const Set<const VariantEvent*> 
		     &coveredEvents);
  void maskLeft(PriorMask &,const Edge &);
  void maskLeftGT(PriorMask &mask,Signal *,const Edge &);
  void maskLeftAG(PriorMask &mask,Signal *,const Edge &);
  void maskLeftATG(PriorMask &mask,Signal *,const Edge &);
  void maskRight(PriorMask &,const Edge &);
  void maskRightGT(PriorMask &mask,Signal *,const Edge &);
  void maskRightAG(PriorMask &mask,Signal *,const Edge &);
  void maskEvents(PriorMask &,const Edge &,
		  const Set<const VariantEvent*> &coveredEvents);
  void regionOverlapping(int pos,ContentType &,Interval &);
  void mask(const Interval &,PriorMask &);
  void computePrior(Edge &,const PriorMask &mask);
  float computePrior(const Labeling &altLabels,int offset,const PriorMask &mask);
  void initLabeling(const Edge &,Labeling &);
  void initExonLabeling(int startingPhase,Labeling &);
};


#endif
