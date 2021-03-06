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
#endif
#include "BOOM/VectorSorter.H"
#include "BOOM/ListQueue.H"
#include "SignalStreamBuilder.H"

CIA::CIA(const BOOM::String &PROGRAM_NAME,
	 const BOOM::String &VERSION,EdgeFactory &edgeFactory,
	 int &transcriptId,VariantEvents &events,const String &
	 projectedGFF,const String &labelFile)
  : GeneZilla(PROGRAM_NAME,VERSION,edgeFactory,transcriptId),
    signalLabelingProfiles(NumContentTypes), priorWeight(0), 
    events(events), refAnno(NULL), projectedGFF(projectedGFF),
    labelFile(labelFile), constraints(NULL)
{
  // ctor
}



CIA::~CIA()
{
  delete refAnno;
}



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
  else {
    BOOM::Vector<IsochoreInterval> &intervals=
      isochoreIntervals[substrateId];
    if(intervals.size()>0) {
      IsochoreInterval &interval=intervals[0];
      gcContent=interval.GC;
      nextIsochoreIndex=1;
      if(intervals.size()>1) nextIsochoreInterval=intervals[1];
      else nextIsochoreInterval.begin=-1;
    }
    else gcContent=getGCcontent(substrateString);
  }

  // This is the first chunk we are seeing...so do some initialization
  // first:
  processIsochoreFile(isoFilename,gcContent);
  const String matrixFile=isochore->configFile.lookupOrDie("label-matrix");
  priorWeight=isochore->configFile.getFloatOrDie("prior-weight");
  int maxIntronScan=isochore->configFile.getIntOrDie("max-intron-scan-length");
  int minExonLength=isochore->configFile.getIntOrDie("min-variant-exon-length");
  int minIntronLen=isochore->configFile.getIntOrDie("min-variant-intron-length");
  bool allowSignalGains=isochore->configFile.getBoolOrDie("allow-signal-gains");
  bool allowGainExonBrokenStop=isochore->configFile.getBoolOrDie("allow-gain-exon-broken-stop");
  shouldReweight=!isochore->configFile.getBoolOrDie("no-prior");

  // Load reference annotation
  refAnno=new ReferenceAnnotation(projectedGFF,labelFile,matrixFile,*isochore,
				  substrateString,substrate);
  initSignalLabelingProfiles();

  // Populate the signal stream
  constraints=new ConstraintIntervals(substrateString.length());
  SignalStreamBuilder ssb(*refAnno,events,signalStream,*constraints,newSignals,
			  isochore,maxIntronScan,minExonLength,minIntronLen,
			  allowSignalGains,allowGainExonBrokenStop);

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

  buildParseGraph(seq,str);
  if(shouldReweight) {
    cerr<<"reweighting graph"<<endl;
    reweightGraph();
  }

  BOOM::Stack<SignalPtr> *path=NULL;
  /*
  double parseScore;
  path=parseGraph.findOptimalPath(parseScore);
  generateGff(path,seqLen,parseScore);
  */
  if(dumpGraph) {
    parseGraph.setVertexIndices();
    osGraph<<parseGraph<<endl;
  }

#ifdef REPORT_MEMORY_USAGE
    MemoryProfiler::report("CIA TOTAL MEMORY USAGE: ",cerr);
#endif

  return path;
}



/*
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
      intronScore=score+priorWeight*
	log(refAnno->getMatrix()(refAnno->getLabeling()[pos],LABEL_INTRON));
    }
    for(int phase=0 ; phase<3 ; ++phase) { // Score against the prior labeling
      GeneModelLabel predictedLabel;
      if(isIntron(contentType)) predictedLabel=LABEL_INTRON; // NEVER HAPPENS
      else if(isIntergenic(contentType) || isUTR(contentType)) 
	predictedLabel=LABEL_INTERGENIC;
      else if(isCoding) predictedLabel=getExonLabel(phase);
      else INTERNAL_ERROR;
      float prior=priorWeight*
	log(refAnno->getMatrix()(refAnno->getLabeling()[pos],predictedLabel));
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
*/


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
    while(1) {
      SignalPtr signal=signalStream.detect(pos);
      if(signal) {
	int begin=signal->getConsensusPosition();
	int end=signal->posOfBaseFollowingConsensus();
	bool supported=false;
	if(evidenceFilter)
	  switch(signal->getSignalType()) {
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
	  //scoreSignalPrior(signal); // ### BEING DONE ELSEWHERE
	  linkBack(str,signal);
	  enforceConstraints(signal);
	  enqueue(signal);
	}	
      }
      else break;
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



void CIA::scoreSignalPrior(SignalPtr s)
{
  if(s->getConsensusPosition()<0 || s->getConsensusPosition()>=seqLen) return;
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
	  if(pos>=seqLen) INTERNAL_ERROR; // ###
	  GeneModelLabel priorLabel=refAnno->getLabeling()[pos];
	  float penalty=priorWeight*
	    refAnno->getMatrix()(priorLabel,predictedLabel);
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
  }
}



void CIA::purgeQueues()
{
  for(Vector<ContentSensor*>::iterator cur=isochore->contentSensors.begin(),
        end=isochore->contentSensors.end() ; cur!=end ; ++cur) {
    ContentSensor &contentSensor=**cur;
    BOOM::Set<SignalQueue*> &queues=contentSensor.getSignalQueues();
    for(Set<SignalQueue*>::iterator cur=queues.begin(), end=queues.end();
        cur!= end ; ++cur) {
      SignalQueue &queue=**cur;
      if(queue.getContentType()==INTERGENIC) continue;
      queue.resetQueue(isochore);
    }
  }
}



void CIA::enforceConstraints(SignalPtr signal)
{
  const int pos=signal->getConsensusPosition();
  if(constraints->isConstrained(pos)) purgeQueues();
}



void CIA::reweightGraph()
{
  Set<Signal*> seen;
  ListQueue<Signal*> Q;
  for(Vector<SignalPtr>::iterator cur=leftTermini.begin(), end=
	leftTermini.end() ; cur!=end ; ++cur) Q.enqueue(*cur);
  while(!Q.isEmpty()) {
    Signal *signal=Q.dequeue();
    seen+=signal;
    Set<Edge*> &edges=signal->getEdgesOut();
    for(Set<Edge*>::iterator cur=edges.begin(), end=edges.end() ;
	cur!=end ; ++cur) {
      Edge &edge=**cur;
      reweight(edge);
      Signal *right=edge.getRight();
      if(!seen.isMember(right)) Q.enqueue(right);
    }
  }
}



void CIA::makePriorMask(PriorMask &mask,const Edge &edge,
			bool leftIsNew,bool rightIsNew,
			const Set<const VariantEvent*> &coveredEvents)
{
  if(leftIsNew) maskLeft(mask,edge);
  if(rightIsNew) maskRight(mask,edge);
  if(coveredEvents.size()>0) maskEvents(mask,edge,coveredEvents);
}



void CIA::regionOverlapping(int pos,ContentType &contentType,Interval &interval)
{
  const ContentRegion *region=refAnno->getRegions().regionOverlapping(pos);
  contentType=region->getType();
  interval=region->getInterval();
}



void CIA::mask(const Interval &I,PriorMask &mask)
{
  const int begin=I.getBegin(), end=I.getEnd();
  for(int i=begin ; i<end ; ++i) mask.mask(i);
}



void CIA::maskLeftGT(PriorMask &pmask,Signal *signal,const Edge &edge)
{
  const int consensusPos=signal->getConsensusPosition();
  ContentType regionType; Interval regionInterval;
  regionOverlapping(consensusPos,regionType,regionInterval);
  Interval edgeInterval=edge.getFeatureInterval().intersect(Interval(0,seqLen));
  if(isCoding(regionType)) {
    Interval maskInterval=regionInterval.intersect(edgeInterval);
    mask(maskInterval,pmask);
  }
}



void CIA::maskLeftAG(PriorMask &pmask,Signal *signal,const Edge &edge)
{
  const int consensusPos=signal->getConsensusPosition();
  ContentType regionType; Interval regionInterval;
  regionOverlapping(consensusPos,regionType,regionInterval);
  const Interval edgeInterval=edge.getFeatureInterval().intersect(Interval(0,seqLen));
  if(isIntron(regionType)) {
    Interval maskInterval=regionInterval.intersect(edgeInterval);
    mask(maskInterval,pmask);
  }
}



void CIA::maskLeftATG(PriorMask &pmask,Signal *signal,const Edge &edge)
{
  const int consensusPos=signal->getConsensusPosition();
  ContentType regionType; Interval regionInterval;
  regionOverlapping(consensusPos,regionType,regionInterval);
  const Interval edgeInterval=edge.getFeatureInterval().intersect(Interval(0,seqLen));
  if(isUTR5(regionType)) {
    Interval maskInterval=regionInterval.intersect(edgeInterval);
    mask(maskInterval,pmask);
  }
}



void CIA::maskLeft(PriorMask &mask,const Edge &edge)
{
  const Signal *signal=edge.getLeft();
  const SignalType signalType=signal->getSignalType();
  switch(signalType) {
  case GT:  maskLeftGT(mask,signal,edge); break;
  case AG:  maskLeftAG(mask,signal,edge); break;
  case ATG: maskLeftATG(mask,signal,edge); break;
  }
}



void CIA::maskRight(PriorMask &mask,const Edge &edge)
{
  const Signal *signal=edge.getRight();
  const SignalType signalType=signal->getSignalType();
  switch(signalType) {
  case GT:  maskRightGT(mask,signal,edge); break;
  case AG:  maskRightAG(mask,signal,edge); break;
  }
}



void CIA::maskRightGT(PriorMask &pmask,Signal *signal,const Edge &edge)
{
  const int consensusPos=signal->getConsensusPosition();
  ContentType regionType; Interval regionInterval;
  regionOverlapping(consensusPos,regionType,regionInterval);
  const Interval edgeInterval=edge.getFeatureInterval().intersect(Interval(0,seqLen));
  if(isIntron(regionType)) {
    Interval maskInterval=regionInterval.intersect(edgeInterval);
    mask(maskInterval,pmask);
  }
}



void CIA::maskRightAG(PriorMask &pmask,Signal *signal,const Edge &edge)
{
  const int consensusPos=signal->getConsensusPosition();
  ContentType regionType; Interval regionInterval;
  regionOverlapping(consensusPos,regionType,regionInterval);
  const Interval edgeInterval=edge.getFeatureInterval().intersect(Interval(0,seqLen));
  if(isCoding(regionType)) {
    Interval maskInterval=regionInterval.intersect(edgeInterval);
    mask(maskInterval,pmask);
  }
}



void CIA::maskEvents(PriorMask &pmask,const Edge &edge,
		     const Set<const VariantEvent*> &coveredEvents)
{
  const Interval edgeInterval=edge.getFeatureInterval().intersect(Interval(0,seqLen));
  for(Set<const VariantEvent*>::const_iterator cur=coveredEvents.begin(),
	end=coveredEvents.end() ; cur!=end ; ++cur) {
    const VariantEvent *event=*cur;
    const VariantEventType eventType=event->getEventType();
    const VariantSignalType signalType=event->getSignalType();
    if(eventType==VAR_EVENT_LOSS && signalType==VAR_SIG_TAG && // loss of stop
       edge.getContentType()==INTERNAL_EXON) {
      ContentType regionType; Interval regionInterval;
      const int pos=(event->begin()+event->end())/2;
      regionOverlapping(pos,regionType,regionInterval);
      if(regionType==INTRON) 
	mask(regionInterval.intersect(edgeInterval),pmask);
    }
  }
}



void CIA::reweight(Edge &edge)
{
  // First, establish prior mask
  Interval featureInterval=edge.getFeatureInterval().intersect(Interval(0,seqLen));
  PriorMask mask(featureInterval);
  const bool leftIsNew=newSignals.isMember(edge.getLeft());
  const bool rightIsNew=newSignals.isMember(edge.getRight());
  Set<const VariantEvent*> coveredEvents;
  events.eventsInInterval(featureInterval,coveredEvents);
  makePriorMask(mask,edge,leftIsNew,rightIsNew,coveredEvents);

  // Apply prior on all unmasked regions
  computePrior(edge,mask);
}



float CIA::computePrior(const Labeling &proposedLabels,int offset,const PriorMask &mask)
{
  const Labeling &priorLabels=refAnno->getLabeling();
  const LabelMatrix &M=refAnno->getMatrix();
  float prior=0;
  const int edgeLen=proposedLabels.length();
  for(int altPos=0, refPos=offset ; altPos<edgeLen ; ++altPos, ++refPos) {
    if(refPos>=seqLen) INTERNAL_ERROR;//###
    if(altPos>=seqLen) INTERNAL_ERROR;//###
    //    cout<<altPos<<" "<<refPos<<" "<<seqLen<<endl;
    //    cout<<priorLabels.length()<<" "<<proposedLabels.length()<<endl;
    if(!mask.lookup(refPos)) {
      prior+=M(priorLabels[refPos],proposedLabels[altPos]);
    }
  }
  return prior;
}



void CIA::initExonLabeling(int startingPhase,Labeling &labeling)
{
  int L=labeling.length();
  if(L>seqLen) INTERNAL_ERROR;//###
  int phase=startingPhase;
  for(int i=0 ; i<L ; ++i) {
    labeling[i]=getExonLabel(phase);
    phase=(phase+1)%3;
  }
}



void CIA::initLabeling(const Edge &edge,Labeling &labeling)
{
  const ContentType contentType=edge.getContentType();
  const Interval featureInterval=edge.getFeatureInterval().intersect(Interval(0,seqLen));
  
  switch(contentType) {
  case INITIAL_EXON: 
    initExonLabeling(0,labeling); 
    break;
  case INTERNAL_EXON: 
    INTERNAL_ERROR;
  case FINAL_EXON: 
    initExonLabeling(posmod(-featureInterval.length()),labeling); 
    break;
  case SINGLE_EXON: 
    initExonLabeling(0,labeling); 
    break;
  case INTRON: 
  case UTR5_INTRON:
  case NEG_UTR5_INTRON:
  case UTR3_INTRON:
  case NEG_UTR3_INTRON:
    labeling.setAllTo(LABEL_INTRON); 
    break;
  case INTERGENIC: 
    labeling.setAllTo(LABEL_INTERGENIC); 
    break;
  case UTR5_INITIAL:
  case UTR5_INTERNAL:
  case UTR5_FINAL:
  case UTR5_SINGLE:
  case NEG_UTR5_INITIAL:
  case NEG_UTR5_INTERNAL:
  case NEG_UTR5_FINAL:
  case NEG_UTR5_SINGLE:
    labeling.setAllTo(LABEL_UTR);
    break;
  case UTR3_INITIAL:
  case UTR3_INTERNAL:
  case UTR3_FINAL:
  case UTR3_SINGLE:
  case NEG_UTR3_INITIAL:
  case NEG_UTR3_INTERNAL:
  case NEG_UTR3_FINAL:
  case NEG_UTR3_SINGLE:
    labeling.setAllTo(LABEL_UTR);
    break;
  default: INTERNAL_ERROR; break;
  }
}



void CIA::computePrior(Edge &edge,const PriorMask &mask)
{
  // Compute prior for edge
  const Interval featureInterval=edge.getFeatureInterval().intersect(Interval(0,seqLen));
  const ContentType contentType=edge.getContentType();
  if(edge.getStrand()==REVERSE_STRAND) INTERNAL_ERROR;// ### DEBUGGING
  if(isCoding(contentType)) { // has three phases
    PhasedEdge &phasedEdge=dynamic_cast<PhasedEdge&>(edge);
    for(int phase=0 ; phase<3 ; ++phase) {
      Labeling labeling(featureInterval.length());
      initExonLabeling(0,labeling);
      const float prior=computePrior(labeling,featureInterval.getBegin(),mask);
      const float newScore=phasedEdge.getEdgeScore(phase)+priorWeight*prior;
      //      cout<<"XXX "<<phase<<" "<<phasedEdge.getEdgeScore(phase)<<" => "<<newScore<<endl;
      phasedEdge.setEdgeScore(phase,newScore);      
    }
  }
  else if(contentType==INTRON) { // has three phases
    PhasedEdge &phasedEdge=dynamic_cast<PhasedEdge&>(edge);
    Labeling labeling(featureInterval.length());
    initLabeling(edge,labeling);
    const float prior=computePrior(labeling,featureInterval.getBegin(),mask);
    for(int phase=0 ; phase<3 ; ++phase) {
      const float newScore=phasedEdge.getEdgeScore(phase)+priorWeight*prior;
      //      cout<<"YYY "<<phase<<" "<<phasedEdge.getEdgeScore(phase)<<" => "<<newScore<<endl;
      phasedEdge.setEdgeScore(phase,newScore);      
    }
  }
  else { // phase can be ignored
    NonPhasedEdge &unphasedEdge=dynamic_cast<NonPhasedEdge&>(edge);
    Labeling labeling(featureInterval.length());
    initLabeling(edge,labeling);
    const float prior=computePrior(labeling,featureInterval.getBegin(),mask);
    const float newScore=unphasedEdge.getEdgeScore()+priorWeight*prior;
    //    cout<<"ZZZ "<<unphasedEdge.getEdgeScore()<<" => "<<newScore<<endl;
    unphasedEdge.setEdgeScore(newScore);
  }

  // Compute prior for right vertex
  if(edge.getRight()->getContextWindowEnd()<seqLen)
    scoreSignalPrior(edge.getRight());
}




