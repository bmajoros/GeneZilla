/****************************************************************
 GeneZilla-CRF
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_CRF_H
#include "crf.H"
#include "GZilla.H"
#include "SignalLabelingProfile.H"

/****************************************************************
                            class CRF
 ****************************************************************/

class CRF : public GeneZilla {
public:
  CRF(const BOOM::String &PROGRAM_NAME,const BOOM::String &VERSION,
      EdgeFactory &,int &transcriptId,LabelMatrix &,float priorWeight);
  virtual ~CRF();
  virtual int main(int argc,char *argv[]);
  virtual BOOM::Stack<SignalPtr> * processChunk(const Sequence &seq,
				      const BOOM::String &seqStr,
				      const BOOM::String &isoFilename,
				      const BOOM::String &substrateId,
				      ostream &osGraph,
				      bool dumpGraph,
				      String psaFilename,
				      Labeling &);

protected:
  Labeling priorLabels;
  LabelMatrix labelMatrix;
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
};


#endif
