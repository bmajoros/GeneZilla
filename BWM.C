/****************************************************************
 BWM.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "BWM.H"
#include "BOOM/Array1D.H"
#include "BOOM/BinomialDistribution.H"
#include <iostream>
#include <fstream>
#include <math.h>
#include "ContentSensor.H"

const double log_of_2=log(2.0);
inline float ln(float p) {return log(p)/log_of_2;}



BWM::BWM(GarbageCollector &gc,BOOM::Vector<TrainingSequence*> &sequences,
	 SignalType signalType,int consensusOffset,int consensusLength,
	 float gcContent,float alpha)
  : WMM(gc)
{
  /*
    This constructor performs training of the BWM
   */

  // Set some things in the base class
  setStrand(FORWARD_STRAND);
  setSignalType(signalType);

  // Compute background single-nucleotide probabilities
  float atContent=1-gcContent;
  float AT=atContent/2, GC=gcContent/2;
  BOOM::Array1D<float> background(alphabet.getNumElements());
  background.setAllTo(0);
  Symbol A=alphabet.lookup('A'), T=alphabet.lookup('T');
  Symbol C=alphabet.lookup('C'), G=alphabet.lookup('G');
  background[A]=AT;
  background[T]=AT;
  background[C]=GC;
  background[G]=GC;

  // Allocate the underlying WMM matrix
  float pseudocount=0;
  int n=sequences.size();
  int len=sequences[0]->getLength();
  setSizes(consensusLength,consensusOffset,len);
  int nAlpha=alphabet.getNumElements();
  matrix.resize(len,nAlpha);
  matrix.setAllTo(pseudocount);

  // Iterate through the training sequences & collect counts
  BOOM::FloatArray1D effectiveSize(len);
  effectiveSize.setAllTo(0);
  for(int i=0 ; i<n ; ++i)
    {
      TrainingSequence &seq=*sequences[i];
      int l=seq.getLength();
      if(l!=len) throw BOOM::String("length mismatch in BWM: ")+len+" vs "+l;
      for(int pos=0 ; pos<len ; ++pos)
	{
	  Symbol s=seq[pos];
	  int count=seq.getBoostCount();
	  matrix[pos][s]+=count;
	  effectiveSize[pos]+=count;
	}
    }

  // Perform a binomial test at each position to see if frequencies
  // are significantly different from the background frequencies
  numSignificantPositions=0;
  for(int pos=0 ; pos<len ; ++pos)
    {
      // First, we have to identify the most extreme count
      int mostExtremeDiff=0;
      Symbol mostExtremeSymbol=0;
      int sampleSize=int(effectiveSize[pos]+5/9.0);
      for(Symbol s=0 ; s<nAlpha ; ++s)
	{
	  int expected=int(background[s]*sampleSize+5/9.0);
	  int observed=(int)matrix[pos][s];
	  int diff=abs(expected-observed);
	  if(diff>mostExtremeDiff)
	    {
	      mostExtremeDiff=diff;
	      mostExtremeSymbol=s;
	    }
	}

      // Apply binomial distribution to get P-value for this count
      float p=background[mostExtremeSymbol];
      int expected=int(p*sampleSize+5/9.0);
      int observed=expected+mostExtremeDiff; // +/- are symmetric, so use +
      double P=
	BinomialDistribution::rightTailedPValue(observed,sampleSize,p);

      // If not significantly different from background, then the observed
      // frequencies are probably just a statistical fluctation due to
      // small sample size, so use the background probabilities instead
      if(P>alpha)
	{
	  effectiveSize[pos]=0;
	  for(Symbol s=0 ; s<nAlpha ; ++s)
	    {
	      matrix[pos][s]=background[s]*sampleSize;
	      effectiveSize[pos]+=matrix[pos][s];
	    }
	}
      else 
	{
	  cout<<pos<<" : using OBSERVED frequencies, p="<<P<<endl;
	  ++numSignificantPositions;
	}
    }

  // Normalize counts into probabilities
  for(int pos=0 ; pos<len ; ++pos)
    for(Symbol s=0 ; s<nAlpha ; ++s)
      matrix[pos][s]/=effectiveSize[pos];

  // Convert probabilities to log-probabilities
  convertToLogs();  
}
