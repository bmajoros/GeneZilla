/****************************************************************
 GeneZilla-CIA
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <fstream>
#include "CIA.H"
#include "BOOM/FastaReader.H"
#ifdef REPORT_PROGRESS
#include "BOOM/Progress.H"
#include "BOOM/VectorSorter.H"
#endif

CIA::CIA(const BOOM::String &PROGRAM_NAME,
	 const BOOM::String &VERSION,EdgeFactory &edgeFactory,
	 int &transcriptId,VariantEvents &events,const String &
	 projectedGFF,const String &labelFile)
  : GeneZilla(PROGRAM_NAME,VERSION,edgeFactory,transcriptId),
    signalLabelingProfiles(NumContentTypes), priorWeight(0), 
    events(events), refAnno(NULL), projectedGFF(projectedGFF),
    labelFile(labelFile)
{
  // ctor
}



CIA::~CIA()
{
  delete refAnno;
}



/*
int CIA::main(int argc,char *argv[])
  {
    // Process command line
    BOOM::CommandLine cmd(argc,argv,"r:");
    if(cmd.numArgs()!=2)
      throw string(
"\ngenezilla <*.iso> <*.fasta> [-r <*.gff>]\n\
          where -r <*.gff> specifies a GFF file to load and score\n\n");
    BOOM::String isochoreFilename=cmd.arg(0);
    BOOM::String fastaFilename=cmd.arg(1);

    cerr << "Loading sequence from " << fastaFilename 
	 << "........................" << endl;
    alphabet=DnaAlphabet::global();

    BOOM::FastaReader fastaReader(fastaFilename);
    BOOM::String defline, substrateString;
    fastaReader.nextSequence(defline,substrateString);
    Sequence substrate(substrateString,DnaAlphabet::global());

    seqLen=substrate.getLength();

//    if(cmd.option('r')) 
//      loadGFF(cmd.optParm('r'),substrate,substrateString);

    cerr << "Running the CIA algorithm..." << endl;
    cerr << "\tsequence length=" << substrate.getLength() << endl;
    ofstream osGraph;
    BOOM::Stack<SignalPtr> *path=mainAlgorithm(substrate,substrateString,
					       osGraph,false,"");
    delete path;

#ifdef REPORT_MEMORY_USAGE
    MemoryProfiler::report("TOTAL MEMORY USAGE: ",cerr);
#endif

    return 0;
  }
*/


BOOM::Stack<SignalPtr> * CIA::processChunk(const Sequence &substrate,
					   const BOOM::String &substrateString,
					   const BOOM::String &isoFilename,
					   const BOOM::String &substrateId,
					   ostream &osGraph,
					   bool dumpGraph,
					   String psaFilename)
{
  seq=&substrate;
  seqStr=&substrateString;
  seqLen=substrate.getLength();
  this->substrateId=substrateId;
  
  nextIsochoreInterval.begin=-1;
  if(!isochoreIntervals.isDefined(substrateId))
    gcContent=getGCcontent(substrateString);
  else
    {
      BOOM::Vector<IsochoreInterval> &intervals=
	isochoreIntervals[substrateId];
      if(intervals.size()>0)
	{
	  IsochoreInterval &interval=intervals[0];
	  gcContent=interval.GC;
	  nextIsochoreIndex=1;
	  if(intervals.size()>1)
	    nextIsochoreInterval=intervals[1];
	  else
	    nextIsochoreInterval.begin=-1;
	}
      else
	gcContent=getGCcontent(substrateString);
    }

  // This is the first chunk we are seeing...so do some initialization
  // first:
  cerr << "Processing config file..." << endl;
  processIsochoreFile(isoFilename,gcContent);
  initSignalLabelingProfiles();
  const String matrixFile=isochore->configFile.lookupOrDie("label-matrix");
  priorWeight=isochore->configFile.getFloatOrDie("prior-weight");

  // Load reference annotation
  refAnno=new ReferenceAnnotation(projectedGFF,labelFile,matrixFile,*isochore,
				  substrateString,substrate);

  cerr << "\tsequence length=" << seqLen << endl;
  return mainAlgorithm(substrate,substrateString,osGraph,dumpGraph,
		       psaFilename);
}


BOOM::Stack<SignalPtr> * CIA::mainAlgorithm(const Sequence &seq,
						  const BOOM::String &str,
						  ostream &osGraph,
						  bool dumpGraph,
						  String psaFilename)
{
  // Compute cumulative intergenic score at each base
  const char *charPtr=str.c_str();
  intergenicSums.resize(seqLen);
  computeIntergenicSums(seq,str,charPtr);
  if(psaFilename.length()>0) {
    ofstream os(psaFilename.c_str());
    for(int i=0 ; i<seqLen ; ++i)
      os<<intergenicSums[i]<<endl;
  }

#ifndef EXPLICIT_GRAPHS
  // Instantiate one signal of each type at the left terminus to act as
  // anchors to which real signals can link back
  instantiateLeftTermini();

#ifdef REPORT_PROGRESS
  BOOM::Progress progress;
  progress.start(seqLen);
#endif

  // Make a single left-to-right pass across the sequence
  for(int pos=0 ; pos<seqLen ; ++pos, ++charPtr)
    {
      if(pos==nextIsochoreInterval.begin) crossIsochoreBoundary(pos);
      Symbol base=seq[pos];
      
      // Check whether any signals occur here
      BOOM::Vector<SignalSensor*>::iterator cur=
	isochore->signalSensors.begin(), end=isochore->signalSensors.end();
      for(; cur!=end ; ++cur )
	{
	  SignalSensor &sensor=**cur;
	  if(pos+sensor.getContextWindowLength()>seqLen) continue;
#ifdef FORCE_SPECIFIC_SIGNALS
	  SignalPtr signal=
	    (forcedSignalCoords.isMember(pos+sensor.getConsensusOffset()) ? 
	     sensor.detectWithNoCutoff(seq,str,pos) :
	     sensor.detect(seq,str,pos));
#else
	  SignalPtr signal=sensor.detect(seq,str,pos);
#endif
	  if(signal)
	    {
	      scoreSignalPrior(signal);

	      // Find optimal predecessor for this signal in all 3 phases
	      linkBack(str,signal);

	      // Add this signal to the appropriate queue(s)
	      enqueue(signal);
	    }
	}

      // Check for stop codons & terminate reading frames when they 
      // are found.  This check lags behind by two bases so that any 
      // stop codon we find won't overlap with a signal yet to be 
      // identified (consider TAGT=TAG+GT; the TAG should not stop 
      // any reading frames for the GT because when GT is used as a 
      // donor only the TA would be present during translation)
      if(pos>1) handleStopCodons(str,pos-2);

      // Propagate scores of all non-eclipsed signals up to this base
      updateAccumulators(seq,str,pos,base,*charPtr);

#ifdef REPORT_PROGRESS
      if((pos+1)%250000==0) cerr<<progress.getProgress(pos)<<endl;
#endif
    }

  // Instantiate an anchor signal of each type at the right terminus
  // and link them back, to complete the dynamic programming evaluation:
  double parseScore;
  BOOM::Stack<SignalPtr> *path=
    instantiateRightTermini(str,seqLen,parseScore);

  // Output gene prediction in GFF format
  generateGff(path,seqLen,parseScore);
#endif

#ifdef EXPLICIT_GRAPHS
  //###DEBUGGING: do the prediction using the graph instead of the "trellis"
  //const char *charPtr=str.c_str();
  //intergenicSums.resize(seqLen);
  //computeIntergenicSums(seq,str,charPtr);
  buildParseGraph(seq,str);
  double parseScore;
  BOOM::Stack<SignalPtr> *path=parseGraph.findOptimalPath(parseScore);
  generateGff(path,seqLen,parseScore);
  if(dumpGraph) {
    parseGraph.setVertexIndices();
    osGraph<<parseGraph<<endl;
  }
#endif

#ifdef REPORT_MEMORY_USAGE
    MemoryProfiler::report("CIA TOTAL MEMORY USAGE: ",cerr);
#endif

  return path;
}



void CIA::updateAccumulators(const Sequence &seq,
				  const BOOM::String &str,
				  int pos,Symbol base,char c)
{
  double score, intronScore, scorePhase0, scorePhase1, scorePhase2;
  for(Vector<ContentSensor*>::iterator cur=isochore->contentSensors.begin(),
	end=isochore->contentSensors.end() ; cur!=end ; ++cur) {
    ContentSensor &contentSensor=**cur;
    ContentType contentType=contentSensor.getContentType();
    bool isCoding=contentSensor.isCoding();
    if(isCoding)
      contentSensor.scoreSingleBase(seq,str,pos,base,c,scorePhase0,
				    scorePhase1,scorePhase2);
    else {
      score=contentSensor.scoreSingleBase(seq,str,pos,base,c);
      intronScore=
	score+priorWeight*log(labelMatrix(refAnno->getLabeling()[pos],LABEL_INTRON));
    }
    for(int phase=0 ; phase<3 ; ++phase) { // Score against the prior labeling
      GeneModelLabel predictedLabel;
      if(isIntron(contentType)) predictedLabel=LABEL_INTRON; // NEVER HAPPENS
      else if(isIntergenic(contentType) || isUTR(contentType)) 
	predictedLabel=LABEL_INTERGENIC;
      else if(isCoding) predictedLabel=getExonLabel(phase);
      else INTERNAL_ERROR;
      float prior=priorWeight*log(labelMatrix(refAnno->getLabeling()[pos],predictedLabel));
      if(predictedLabel==LABEL_INTRON) cout<<"\t"<<score+prior<<endl;
      if(priorWeight==0) prior=0;
      if(isNaN(prior)) INTERNAL_ERROR; // ###
      if(getStrand(contentType)!=PLUS_STRAND) prior=NEGATIVE_INFINITY; // ###
      if(isCoding) switch(phase) {
	case 0: scorePhase0+=prior; break;
	case 1: scorePhase1+=prior; break;
	case 2: scorePhase2+=prior; break;
      }
      else if(phase==0) score+=prior;
    }
    BOOM::Set<SignalQueue*> &queues=contentSensor.getSignalQueues();
    for(BOOM::Set<SignalQueue*>::iterator cur=queues.begin(), end=queues.end();
	cur!= end ; ++cur) {
      SignalQueue &queue=**cur;
      if(isCoding) 
	queue.addToAccumulator(scorePhase0,scorePhase1,scorePhase2,pos);
      else if(isIntron(queue.getContentType()))
	queue.addToAccumulator(intronScore);
      else queue.addToAccumulator(score);
    }
  }
}



#ifdef EXPLICIT_GRAPHS
void CIA::buildParseGraph(const Sequence &seq,const BOOM::String &str)
{
  // Instantiate one signal of each type at the left terminus to act as
  // anchors to which real signals can link back
  instantiateLeftTermini();

  // Make a single left-to-right pass across the sequence
  intergenicSums.resize(seqLen);
  const char *charPtr=str.c_str();
  computeIntergenicSums(seq,str,charPtr);
  for(int pos=0 ; pos<seqLen ; ++pos, ++charPtr) {
    Symbol base=seq[pos];
      
    // Check whether any signals occur here
    BOOM::Vector<SignalSensor*>::iterator cur=
      isochore->signalSensors.begin(),
      end=isochore->signalSensors.end();
    for(; cur!=end ; ++cur ) {
      SignalSensor &sensor=**cur;
      if(pos+sensor.getContextWindowLength()>seqLen) continue;
      
      SignalPtr signal=sensor.detect(seq,str,pos);
#ifdef FORCE_SPECIFIC_SIGNALS
      if(!signal && forcedSignalCoords.isMember
	 (pos+sensor.getConsensusOffset()))
	signal=sensor.detectWithNoCutoff(seq,str,pos);
#endif
      if(signal) {
	int begin=signal->getConsensusPosition();
	int end=signal->posOfBaseFollowingConsensus();
	bool supported=false;
	if(evidenceFilter)
	  switch(signal->getSignalType()) 
	    {
	    case ATG:
	    case TAG:
	    case NEG_ATG:
	    case NEG_TAG:
	      supported=evidenceFilter->codingSignalSupported(begin,end);
	      break;
	    case GT:
	    case NEG_AG:
	      supported=evidenceFilter->spliceOutSupported(begin);
	      break;
	    case AG:
	    case NEG_GT:
	      supported=evidenceFilter->spliceInSupported(end);
	      break;
	    }
	else supported=true;

	if(supported) {
	  scoreSignalPrior(signal);

	  // Find optimal predecessor for this signal in all 3 phases
	  linkBack(str,signal);
		
	  // Add this signal to the appropriate queue(s)
	  enqueue(signal);

	}	
      }
    }
      
    // Check for stop codons & terminate reading frames when they are 
    // found.
    // This check lags behind by two bases so that any stop codon we find 
    // won't overlap with a signal yet to be identified 
    // (consider TAGT=TAG+GT; the TAG should not stop any reading frames 
    // for the GT because when GT is used as a donor only the TA would be 
    // present during translation)
    if(pos>1) handleStopCodons(str,pos-2);

    // Propagate scores of all non-eclipsed signals up to this base
    updateAccumulators(seq,str,pos,base,*charPtr);
  }

  // Instantiate an anchor signal of each type at the right terminus
  // and link them back, to complete the dynamic programming evaluation:
  double parseScore;
  BOOM::Stack<SignalPtr> *path=
    instantiateRightTermini(str,seqLen,parseScore);
  delete path;

  // Run the garbage collector & build the parse graph
  BOOM::Vector<SignalPtr>::iterator rCur=rightTermini.begin(), 
    rEnd=rightTermini.end();
  for(; rCur!=rEnd ; ++rCur) {
    SignalPtr s=*rCur;
    garbageCollector.markLeft(s);
  }
  BOOM::Vector<SignalPtr>::iterator lCur=leftTermini.begin(), 
    lEnd=leftTermini.end();
  for(; lCur!=lEnd ; ++lCur) {
    SignalPtr s=*lCur;
    garbageCollector.markRight(s);
  }
  garbageCollector.sweep();
  BOOM::Set<SignalPtr>::iterator sCur=garbageCollector.signalsBegin(),
    sEnd=garbageCollector.signalsEnd();
  while(sCur!=sEnd) {
    SignalPtr s=*sCur;
    if(!useSignalScores) s->dropSignalScores();
    if(!useContentScores) s->dropContentScores();
    parseGraph.addVertex(s);
    ++sCur;
    garbageCollector.drop(s);
  }
}
#endif



void CIA::scoreSignalPrior(SignalPtr s)
{
  SignalType t=s->getSignalType();
  SignalLabelingProfile &profile=signalLabelingProfiles[t];
  const int wBegin=s->getContextWindowPosition();
  const int L=s->getContextWindowLength();
  const int wEnd=wBegin+L;
  Set<ContentType> &queues=s->belongsInWhichQueues();
  for(Set<ContentType>::iterator cur=queues.begin(), end=queues.end() ;
      cur!=end ; ++cur) {
    ContentType t=*cur;
    Propagator &prop=s->getPropagator(t);
    for(int phase=0 ; phase<3 ; ++phase) {
      float score=0;
      if(priorWeight>0)
	for(int pos=wBegin ; pos<wEnd ; ++pos) {
	  GeneModelLabel predictedLabel=profile.getLabel(phase,pos-wBegin);
	  GeneModelLabel priorLabel=refAnno->getLabeling()[pos];
	  float penalty=priorWeight*log(labelMatrix(priorLabel,predictedLabel));
	  score+=penalty;
	}
      prop[phase]+=score;
    }
  }
}




void CIA::initSignalLabelingProfiles()
{
  BOOM::Vector<SignalSensor*>::iterator cur=
    isochore->signalSensors.begin(), end=isochore->signalSensors.end();
  for(; cur!=end ; ++cur ) {
    SignalSensor &sensor=**cur;
    signalLabelingProfiles[sensor.getSignalType()]=
      SignalLabelingProfile(sensor);
    //cout<<sensor.getSignalType()<<endl;
    //cout<<signalLabelingProfiles[sensor.getSignalType()]<<endl<<endl;
  }
}




