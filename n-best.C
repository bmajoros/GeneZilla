/****************************************************************
 n-best.C : find N best predictions (isoforms)
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "BOOM/Constants.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Array1D.H"
#include "BOOM/Vector.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/FixedSizePriorityQueue.H"
#include "BOOM/Time.H"
#include "BOOM/Stack.H"
#include "BOOM/WigBinary.H"
#include "LightGraph.H"
#include "TrellisLink.H"
using namespace std;
using namespace BOOM;


class Application {
public:
  Application();
  void main(int argc,char *argv[]);
  void buildTrellis(LightGraph &,int N,BOOM::Vector<TrellisLink> &);
protected:
  int autoNmode;
  Array1D< Array1D<TrellisLink> > links;
  int queueCapacity;
  TrellisLinkComparator cmp;
  bool wantIntrons;
  WigBinary *wig;
  Vector<WigInterval> nonzeroRegions;
  Array1D<bool> edgeEmitted;
  int numSupportedIntrons, supportedIntronsEmitted;
  int numSupportedExons, supportedExonsEmitted;
  int transcriptNum;
  int posmod(int x);
  bool autoAssess(Stack<TrellisLink*> &);
  int countSupportedIntrons(LightGraph &);
  int countSupportedExons(LightGraph &);
  void traceback(TrellisLink *endLink,Stack<TrellisLink*> &path);
  void emitHeader();
  bool emit(Stack<TrellisLink*> &,LightGraph &,int parseNum);
  void emitFooter(LightGraph &G,BOOM::Vector<TrellisLink> &parseList);
  void predictUTR(ContentType,Strand,int edgeBegin,int edgeEnd,
		  const String &contigID,const String &parseID,
		  const String &transcriptID,const String &programName,
		  float transcriptScore);
  void emitUTR(String utrType,Strand strand,int utrBegin,
	       int utrEnd,const String &contigID,
	       const String &parseID,const String &transcriptID,
	       const String &programName,float transcriptScore);
};



int main(int argc,char *argv[]) {
  try {	Application app;  app.main(argc,argv); }
  catch(const char *p) { cerr << p << endl; return -1; }
  catch(const string &msg) { cerr << msg.c_str() << endl; return -1; }
  catch(const exception &e) {
    cerr << "STL exception caught in main:\n" << e.what() << endl;
    return -1;
  }
  catch(...) {
    cerr << "Unknown exception caught in main" << endl;
    return -1;
  }
  return 0;
}



Application::Application()
  : wantIntrons(false), queueCapacity(10), wig(NULL), transcriptNum(0)
{
}



void Application::main(int argc,char *argv[])
{
  // Process command line
  BOOM::CommandLine cmd(argc,argv,"q:u:ia:");
  if(cmd.numArgs()!=2)
    throw BOOM::String("\n\
n-best [options] <in.graph> <N> \n\
    -q queue-capacity\n\
    -u pileup-file : infer UTR using RNA-seq reads\n\
    -i : emit introns\n\
    -a mode : automatically determine the number of predictions to make\n\
              mode 0 : disable this feature\n\
              mode 1 : based on introns only\n\
              mode 2 : based on introns and internal exons [DEPRECATED]\n\
              mode 3 : based on introns, but predict more\n\
");
  const String inGraphFile=cmd.arg(0);
  int N=cmd.arg(1).asInt();
  if(cmd.option('q')) queueCapacity=cmd.optParm('q').asInt();
  wantIntrons=cmd.option('i');
  const bool autoN=cmd.option('a');
  //if(autoN && N<100) N=100;
  //if(autoN && N<50) N=50;
  //if(autoN && N<30) N=30;
  if(autoN) queueCapacity=N;
  autoNmode=autoN ? cmd.optParm('a').asInt() : 0;
  if(autoNmode<0 || autoNmode>3) throw "bad argument to -a option";

  // Load graph
  LightGraph *G=new LightGraph(inGraphFile);
  const int numEdges=G->getNumEdges();
  if(autoNmode>0) {
    numSupportedIntrons=countSupportedIntrons(*G);
    supportedIntronsEmitted=0;
  }	
  if(autoNmode==2) {
    numSupportedExons=countSupportedExons(*G);
    supportedExonsEmitted=0;
  }
  edgeEmitted.resize(G->getLargestEdgeID()+1);//numEdges);
  edgeEmitted.setAllTo(false);

  // Load pileup file
  if(cmd.option('u')) {
    wig=new WigBinary(cmd.optParm('u'));
    wig->regionsAbove(0,nonzeroRegions);
  }

  // Build trellis
  BOOM::Vector<TrellisLink> termini;
  buildTrellis(*G,N,termini);

  // Construct paths & emit GFF
  emitHeader();
  int parseNum=1;
  const int numTermini=termini.size();
  for(int i=0 ; i<numTermini ; ++i) {
    Stack<TrellisLink*> path;
    traceback(&termini[i],path);
    bool shouldEmit=autoN ? autoAssess(path) : true;
    if(shouldEmit || autoNmode==3)
      if(emit(path,*G,parseNum)) ++parseNum;
      else break;
    if((autoNmode==1 || autoNmode==3) 
       && supportedIntronsEmitted>=numSupportedIntrons) break;
    else if(autoNmode==2 && supportedIntronsEmitted>=numSupportedIntrons
	    && supportedExonsEmitted>=numSupportedExons) break;
  }
  emitFooter(*G,termini);
}



int Application::countSupportedIntrons(LightGraph &G)
{
  int n=0;
  int numEdges=G.getNumEdges();
  for(int i=0 ; i<numEdges ; ++i) {
    LightEdge *e=G.getEdge(i);
    if(dropStrand(e->getEdgeType())==INTRON && e->isSupported()) ++n;
  }
  return n;
}



int Application::countSupportedExons(LightGraph &G)
{
  int n=0;
  int numEdges=G.getNumEdges();
  for(int i=0 ; i<numEdges ; ++i) {
    LightEdge *e=G.getEdge(i);
    if(isCoding(dropStrand(e->getEdgeType())) && e->isSupported()) ++n;
  }
  return n;
}



bool Application::autoAssess(Stack<TrellisLink*> &path)
{
  bool retval=false;
  for(Stack<TrellisLink*>::iterator cur=path.begin(), end=path.end() ;
      cur!=end ; ++cur) {
    TrellisLink *pl=*cur;
    LightEdge *e=pl->getEdge();
    if(!e) continue;
    const ContentType edgeType=dropStrand(e->getEdgeType());
    if(edgeType==INTRON && e->isSupported() && !edgeEmitted[e->getID()]) {
      retval=true;
      edgeEmitted[e->getID()]=true;
      ++supportedIntronsEmitted;
    }
    else if(isCoding(edgeType) && e->isSupported() && 
	    !edgeEmitted[e->getID()]) {
      retval=true;
      edgeEmitted[e->getID()]=true;
      ++supportedExonsEmitted;
    }
  }
  return retval;
}



int Application::posmod(int x) 
{
  int f=x%3;
  if(f>=0) return f;
  return f+3;
}



void Application::buildTrellis(LightGraph &G,int N,
			       BOOM::Vector<TrellisLink> &termini)
{
  const int numVertices=G.getNumVertices();
  links.resize(numVertices);
  if(numVertices==0) return;
  for(int i=0 ; i<numVertices ; ++i) {
    FixedSizePriorityQueue<TrellisLink> *Q[3];
    for(int j=0 ; j<3 ; ++j) 
      Q[j]=new FixedSizePriorityQueue<TrellisLink>(queueCapacity,cmp);
    LightVertex *signal=G.getVertex(i);
    if(!signal) continue;
    SignalType signalType=dropStrand(signal->getSignalType());
    Strand signalStrand=signal->getStrand();
    int defaultSignalPhase=signalStrand==FORWARD_STRAND ? 0 : 2;
    Vector<LightEdge*> &edges=signal->getEdgesIn();
    int numEdges=edges.size();
    if(numEdges==0) { //left-terminus:
      TrellisLink link(NULL,NULL);
      Q[0]->insert(link);
    }
    for(Vector<LightEdge*>::iterator cur=edges.begin(), end=edges.end() ; 
	cur!=end ; ++cur) {
      LightEdge *currentEdge=*cur;
      LightVertex *pred=currentEdge->getLeft();
      Array1D<TrellisLink> &predLinks=links[pred->getID()];
      const int numPredLinks=predLinks.size();
      if(!currentEdge->isIntergenic()) { //exon or intron
	if(currentEdge->isCoding()) { //exon
	  for(int j=0 ; j<numPredLinks ; ++j) {
	    TrellisLink &predLink=predLinks[j];
	    const int phase=predLink.getPhase();
	    const int signalPhase=currentEdge->propagateForward(phase);
	    if(signalStrand==FORWARD_STRAND) {
	      if((signalType==TAG || signalType==ATG) && signalPhase!=0)
		continue; }
	    else // REVERSE_STRAND
	      if((signalType==TAG || signalType==ATG) && signalPhase!=2)
		continue;
	    double score=currentEdge->getScore(phase)+
	      predLink.getScore();
	    if(isFinite(score)) {
	      TrellisLink link(&predLink,currentEdge,signalPhase,score);
	      Q[signalPhase]->insert(link);
	    }
	  }
	}
	else { //intron
	  for(int j=0 ; j<numPredLinks ; ++j) {
	    TrellisLink &predLink=predLinks[j];
	    const int phase=predLink.getPhase();
	    const double score=currentEdge->getScore(phase)+
	      predLink.getScore();
	    if(isFinite(score)) {
	      TrellisLink link(&predLink,currentEdge,phase,score);
	      Q[phase]->insert(link);
	    }	
	  }
	}
      }
      else { //intergenic or UTR
	for(int j=0 ; j<numPredLinks ; ++j) {
	  TrellisLink &predLink=predLinks[j];
	  const int phase=defaultSignalPhase;
	  const double score=currentEdge->getScore(0)+predLink.getScore();
	  if(isFinite(score)) {
	    TrellisLink link(&predLink,currentEdge,phase,score);
	    Q[phase]->insert(link);
	  }
	}
      }
    }

    // Copy selected links into signal's link set
    Array1D<TrellisLink> &linkSet=links[i];
    int size=0;
    for(int f=0 ; f<3 ; ++f) size+=Q[f]->getNumElements();
    linkSet.resize(size);
    int j=0;
    for(int f=0 ; f<3 ; ++f)
      for(FixedSizePriorityQueue<TrellisLink>::iterator cur=Q[f]->begin(), 
	    end=Q[f]->end() ; cur!=end ; ++cur)
	linkSet[j++]=*cur;
    for(int j=0 ; j<3 ; ++j) delete Q[j]; // this is NOT the problem
  }
  
  // Find all of the right-terminal links
  for(int rCur=numVertices-1; rCur>=0; --rCur) {
    LightVertex *rt=G.getVertex(rCur);
    if(!rt) continue;
    if(!rt->getEdgesOut().isEmpty()) break;
    Array1D<TrellisLink> &rtLinks=links[rCur];
    const int numRtLinks=rtLinks.size();
    for(int j=0 ; j<numRtLinks ; ++j)
      termini.push_back(rtLinks[j]);
  }
  
  // Return links representing the N best parses
  VectorSorter<TrellisLink> sorter(termini,cmp);
  sorter.sortDescendInPlace();
  if(N<termini.size()) termini.resize(N);
}



void Application::traceback(TrellisLink *endLink,Stack<TrellisLink*> &path) 
{
  TrellisLink *currentLink=endLink;
  while(currentLink) {
    path.push(currentLink);
    currentLink=currentLink->getPred();
  }
}



void Application::emitHeader()
{
  cout<<"##gff-version 2"<<endl;
  cout<<"##source-version Decoder 1.0"<<endl;
  cout<<"##date "<<getDateAndTime()<<endl;
  cout<<"##Type DNA"<<endl;
}



bool Application::emit(Stack<TrellisLink*> &pStack,LightGraph &G,int parseNum)
{
  const String contigID=G.getSubstrate();
  const String programName="RSVP";
  cout<<"\n# parse "<<parseNum<<":"<<endl;
  //int transcriptNum=0;
  Vector<TrellisLink*> pList;
  while(!pStack.isEmpty()) {
    TrellisLink *pl=pStack.pop();
    pList.push_back(pl);
  }
  TrellisLink *end=pList[pList.size()-1];
  const float transcriptScore=exp(end->getScore()/G.getSubstrateLength())/.25;
  //if(transcriptScore<1) return false;
  for(Vector<TrellisLink*>::iterator cur=pList.begin(), end=pList.end() ;
      cur!=end ; ++cur) {
    TrellisLink *pl=*cur;
    LightEdge *e=pl->getEdge();
    if(!e) continue;
    const ContentType edgeType=dropStrand(e->getEdgeType());
    if(edgeType==INTERGENIC) continue;
    if(edgeType==INTRON && !wantIntrons) continue;
    String etString=contentTypeNiceString(edgeType);
    int startPos=e->getBegin()+1;
    int endPos=e->getEnd();
    Strand strand=e->getStrand();
    int length=endPos-startPos+1;
    int phase=e->propagateBackward(pl->getPhase());
    int displayPhase=
      strand==FORWARD_STRAND ? phase : posmod(pl->getPhase()+1);
    String scoreString=String(exp(e->getScore(phase)/length)/.25);
    if(((edgeType==INITIAL_EXON || edgeType==SINGLE_EXON) && 
	strand==FORWARD_STRAND)
	|| ((edgeType==FINAL_EXON || edgeType==SINGLE_EXON) && 
	    strand==REVERSE_STRAND))
      ++transcriptNum;
    String transcriptID=String("transcript_id=")+transcriptNum+";";
    String parseString=String("parse=")+parseNum+";";
    int support=e->isSupported() ? 1 : 0;
    cout<<contigID << "\t" << programName << "\t" << etString << 
      "\t" << startPos << "\t" << endPos << "\t" << scoreString << "\t" << 
      strand << "\t" << displayPhase << "\t" << transcriptID << parseString <<
      "transcript_score="<<transcriptScore<<";"<<"sup="<<support<<";"<<endl;
    predictUTR(edgeType,strand,startPos-1,endPos,contigID,parseString,
	       transcriptID,programName,transcriptScore);
  }
  return true;
}



void Application::emitFooter(LightGraph &G,
			     BOOM::Vector<TrellisLink> &parseList) 
{
  cout<<"\n# sequence length: "<<G.getSubstrateLength()<<endl;
}



void Application::emitUTR(String utrType,Strand strand,int utrBegin,
			  int utrEnd,const String &contigID,
			  const String &parseID,const String &transcriptID,
			  const String &programName,float transcriptScore)
{
  float score=0;// ### wig->ave(utrBegin,utrEnd);
  cout<<contigID << "\t" << programName << "\t" << utrType << 
    "\t" << utrBegin+1 << "\t" << utrEnd << "\t" << score << 
    "\t" << strand << "\t.\t" << transcriptID << parseID <<
    "transcript_score="<<transcriptScore<<";"<<endl;
}


void Application::predictUTR(ContentType edgeType,Strand strand,int edgeBegin,
			     int edgeEnd,const String &contigID,
			     const String &parseID,const String &transcriptID,
			     const String &programName,float transcriptScore)
{
  Vector<WigInterval>::iterator cur=nonzeroRegions.begin(),
    end=nonzeroRegions.end();
  for(; cur!=end ; ++cur) {
    WigInterval wi=*cur;
    //cout<<wi.begin<<" "<<wi.end<<endl;
    if(wi.begin<edgeEnd && edgeBegin<wi.end) {
      switch(edgeType)
	{
	case INITIAL_EXON: {
	  if(strand==FORWARD_STRAND) {
	    int utrBegin=wi.begin, utrEnd=edgeBegin;
	    if(utrBegin+10<utrEnd)
	      emitUTR("5'UTR",strand,utrBegin,utrEnd,contigID,parseID,
		      transcriptID,programName,transcriptScore);
	  }
	  else {
	    int utrBegin=edgeEnd, utrEnd=wi.end;
	    if(utrBegin+10<utrEnd)
	      emitUTR("5'UTR",strand,utrBegin,utrEnd,contigID,parseID,
		      transcriptID,programName,transcriptScore);
	  }
	  break;}

	case FINAL_EXON: {
	  if(strand==FORWARD_STRAND) {
	    int utrBegin=edgeEnd, utrEnd=wi.end;
	    if(utrBegin+10<utrEnd)
	      emitUTR("3'UTR",strand,utrBegin,utrEnd,contigID,parseID,
		      transcriptID,programName,transcriptScore);
	  }
	  else {
	    int utrBegin=wi.begin, utrEnd=edgeBegin;
	    if(utrBegin+10<utrEnd)
	      emitUTR("3'UTR",strand,utrBegin,utrEnd,contigID,parseID,
		      transcriptID,programName,transcriptScore);
	  }
	  break;}

	case SINGLE_EXON: {
	  if(strand==FORWARD_STRAND) {
	    int utrBegin=wi.begin, utrEnd=edgeBegin;
	    if(utrBegin+10<utrEnd)
	      emitUTR("5'UTR",strand,utrBegin,utrEnd,contigID,parseID,
		      transcriptID,programName,transcriptScore);
	    utrBegin=edgeEnd; utrEnd=wi.end;
	    if(utrBegin+10<utrEnd)
	      emitUTR("3'UTR",strand,utrBegin,utrEnd,contigID,parseID,
		      transcriptID,programName,transcriptScore);
	  }
	  else {
	    int utrBegin=edgeEnd, utrEnd=wi.end;
	    if(utrBegin+10<utrEnd)
	      emitUTR("5'UTR",strand,utrBegin,utrEnd,contigID,parseID,
		      transcriptID,programName,transcriptScore);
	    utrBegin=wi.begin; utrEnd=edgeBegin;
	    if(utrBegin+10<utrEnd)
	      emitUTR("3'UTR",strand,utrBegin,utrEnd,contigID,parseID,
		      transcriptID,programName,transcriptScore);
	  }
	  break;}

	case NEG_INITIAL_EXON:{
	  int utrBegin=edgeEnd, utrEnd=wi.end;
	  if(utrBegin+10<utrEnd)
	    emitUTR("5'UTR",strand,utrBegin,utrEnd,contigID,parseID,
		    transcriptID,programName,transcriptScore);
	  break;}

	case NEG_FINAL_EXON:{
	  int utrBegin=wi.begin, utrEnd=edgeBegin;
	  if(utrBegin+10<utrEnd)
	    emitUTR("3'UTR",strand,utrBegin,utrEnd,contigID,parseID,
		    transcriptID,programName,transcriptScore);
	  break;}

	case NEG_SINGLE_EXON:{
	  int utrBegin=edgeEnd, utrEnd=wi.end;
	  if(utrBegin+10<utrEnd)
	    emitUTR("5'UTR",strand,utrBegin,utrEnd,contigID,parseID,
		    transcriptID,programName,transcriptScore);
	  utrBegin=wi.begin; utrEnd=edgeBegin;
	  if(utrBegin+10<utrEnd)
	    emitUTR("3'UTR",strand,utrBegin,utrEnd,contigID,parseID,
		    transcriptID,programName,transcriptScore);
	  break;}
	}
    }
  }
}




