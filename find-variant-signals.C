/****************************************************************
 find-variant-signals.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "GZilla.H"
#include "genezilla.H"
#include "EdgeFactory.H"
#include "BOOM/FastaReader.H"
#include "BOOM/Constants.H"
#include "BOOM/CigarString.H"
using namespace std;

#ifdef EXPLICIT_GRAPHS
#error Please edit the file genezilla.H, comment out the definition of EXPLICIT_GRAPHS, issue a "make clean", and recompile this project
#endif

#ifdef FORCE_SPECIFIC_SIGNALS
#error Please edit the file genezilla.H, comment out the definition of FORCE_SPECIFIC_SIGNALS, issue a "make clean", and recompile this project
#endif

static const char *PROGRAM_NAME="find-variant-signals";
static const char *VERSION="1.0";
Alphabet alphabet;
int frame; // ### CAUTION: this is required by older code; to be removed



/***************************************************************************
                            struct VariantRegion
 ***************************************************************************/
struct VariantRegion {
  int refBegin, refEnd; // zero-based, end is not inclusive
  int altBegin, altEnd; // zero-based, end is not inclusive
  VariantRegion(int refBegin,int refEnd,int altBegin,int altEnd)
    : refBegin(refBegin), refEnd(refEnd), altBegin(altBegin), altEnd(altEnd){}
  VariantRegion(){}
};



/***************************************************************************
                            class Application
 ***************************************************************************/
class Application {
public:
  void AppMain(int argc,char *argv[]);
protected:
  String refSubstrate, altSubstrate;
  bool parseVariantLine(const String &line,VariantRegion &);
  void applySensor(SignalSensor &,const VariantRegion &,Sequence &refSeq,
		   Sequence &altSeq,const String &refSeqStr,const String 
		   &altSeqStr,const CigarAlignment &fw,const CigarAlignment &rev);
  void mapWindow(int refPos,int windowLen,const CigarAlignment &,Set<int> &into,
		 int altLen);
  void emit(SignalPtr,int consensusPosition,const String &substrate,
	    const String &status);
  void emitIndels(const CigarAlignment &);
  int mapConsensus(int refPos,const CigarAlignment &);
};



/***************************************************************************
                                    main()
 ***************************************************************************/
int main(int argc,char *argv[]) {
  try {
    Application app;
    app.AppMain(argc,argv);
  }
  catch(const char *p) {
    cerr << p << endl;
    return -1;
  }
  catch(const string &msg) {
    cerr << msg.c_str() << endl;
    return -1;
  }
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



/***************************************************************************
                                  AppMain()
 ***************************************************************************/
void Application::AppMain(int argc,char *argv[])
{
  // Process command line
  BOOM::CommandLine cmd(argc,argv,"s:i:c:o:I:tP:SCDO");
  if(cmd.numArgs()!=5)
    throw BOOM::String(
    "\nfind-variant-signals <variants.txt> <ref.fasta> <alt.fasta> <ref-alt.cigar> <*.iso>\n");
  const String variantFilename=cmd.arg(0);
  const String refFasta=cmd.arg(1);
  const String altFasta=cmd.arg(2);
  const String cigarFile=cmd.arg(3);
  const String isochoreFilename=cmd.arg(4);
  alphabet=DnaAlphabet::global();
  
  // Load sequences & alignment
  EdgeFactory factory;
  int transcriptId=-1;
  GeneZilla genezilla(PROGRAM_NAME,VERSION,factory,transcriptId);
  BOOM::String defline, refSeqStr, altSeqStr, junk;
  FastaReader::load(refFasta,defline,refSeqStr);
  FastaReader::parseDefline(defline,refSubstrate,junk);
  FastaReader::load(altFasta,defline,altSeqStr);
  FastaReader::parseDefline(defline,altSubstrate,junk);
  Sequence refSeq(refSeqStr,alphabet);
  Sequence altSeq(altSeqStr,alphabet);
  const float gc=genezilla.getGCcontent(refSeqStr);
  IsochoreTable isochores(genezilla.getGC());
  isochores.load(isochoreFilename);
  Isochore *isochore=isochores.getIsochore(gc);
  CigarString cigar; cigar.load(cigarFile);
  CigarAlignment &alignment=*cigar.getAlignment();
  CigarAlignment &revAlign=*alignment.invert(altSeqStr.length());

  // Iterate through variants
  File variantFile(variantFilename);
  VariantRegion region;
  while(!variantFile.eof()) {
    String line=variantFile.getline();
    if(!parseVariantLine(line,region)) continue;

    // Iterate through signal sensors
    for(Vector<SignalSensor*>::iterator cur=isochore->signalSensors.begin(), 
	  end=isochore->signalSensors.end() ; cur!=end ; ++cur) {
      SignalSensor *sensor=*cur;
      if(sensor->getStrand()!=PLUS_STRAND) continue;
      applySensor(*sensor,region,refSeq,altSeq,refSeqStr,altSeqStr,
		  alignment,revAlign);
    }
  }

  // Report indels
  emitIndels(revAlign);
}



bool Application::parseVariantLine(const String &line,VariantRegion &region)
{
  line.trimWhitespace();
  Vector<String> &fields=*line.getFields();
  bool ok=false;
  if(fields.size()==6) {
    region.refBegin=fields[1].asInt();
    region.refEnd=fields[2].asInt();
    region.altBegin=fields[4].asInt();
    region.altEnd=fields[5].asInt();
    ok=true;
  }
  delete &fields;
  return ok;
}



void Application::applySensor(SignalSensor &sensor,const VariantRegion &region,
			      Sequence &refSeq,Sequence &altSeq,const String
			      &refSeqStr,const String &altSeqStr,
			      const CigarAlignment &alignment,
			      const CigarAlignment &revAlign)
{
  const float threshold=sensor.getCutoff();
  const int windowLen=sensor.getContextWindowLength();
  const int refLen=refSeq.getLength(), altLen=altSeq.getLength();

  // Find reference signals missing in the alternate
  {
    int begin=region.refBegin-windowLen+1;
    if(begin<0) begin=0;
    int end=region.refEnd-1;
    if(end+windowLen>refLen) end=refLen-windowLen;
    for(int r=begin ; r<=end ; ++r) {
      SignalPtr signal=sensor.detect(refSeq,refSeqStr,r);
      if(signal) {
	Set<int> altPositions;
	mapWindow(r,windowLen,alignment,altPositions,altLen);
	bool found=false;
	for(Set<int>::const_iterator cur=altPositions.begin(), end=
	      altPositions.end() ; cur!=end ; ++cur) {
	  const int altPos=*cur;
	  SignalPtr altSignal=sensor.detect(altSeq,altSeqStr,altPos);
	  if(altSignal) found=true;
	}
	if(!found) {
	  int altPos=mapConsensus(signal->getConsensusPosition(),alignment);
	  emit(signal,altPos,altSubstrate,"broken");
	}
      }
    }
  }

  // Find alternate signals missing in the reference
  {
    int begin=region.altBegin-windowLen+1;
    if(begin<0) begin=0;
    int end=region.altEnd-1;
    if(end+windowLen>altLen) end=altLen-windowLen;
    for(int r=begin ; r<=end ; ++r) {
      SignalPtr signal=sensor.detect(altSeq,altSeqStr,r);
      if(signal) {
	Set<int> refPositions;
	mapWindow(r,windowLen,revAlign,refPositions,refLen);
	bool found=false;
	for(Set<int>::const_iterator cur=refPositions.begin(), end=
	      refPositions.end() ; cur!=end ; ++cur) {
	  const int refPos=*cur;
	  SignalPtr refSignal=sensor.detect(refSeq,refSeqStr,refPos);
	  if(refSignal) found=true;
	}
	if(!found) emit(signal,signal->getConsensusPosition(),altSubstrate,"new");
      }
    }
  }
}



void Application::mapWindow(int refBegin,int windowLen,
			    const CigarAlignment &alignment,
			    Set<int> &into,int altLen)
{
  /* In order to account for possible indels in the window, we
     map each window position separately, then figure out where
     that would place the beginning of the window.  This gives a
     set of possible starting positions for the signal window in
     the alt sequence, accounting for the fact that indels may
     allow multiple reasonable starting positions */
  for(int i=0 ; i<windowLen ; ++i) {
    const int refPos=refBegin+i;
    const int altPos=alignment[refPos];
    if(altPos==CIGAR_UNDEFINED) continue;
    const int altBegin=altPos-i;
    if(altBegin<0) continue;
    if(altBegin+windowLen>altLen) continue;
    into+=altBegin;
  }
}



void Application::emit(SignalPtr signal,int consPos,const String &substrate,
		       const String &status)
{
  cout<<substrate<<"\t"<<status<<"\t"<<signalTypeToName(signal->getSignalType())
      <<"\t"<<consPos+1<<"\t"<<consPos+signal->getConsensusLength()<<"\t"<<".\t"
      <<signal->getStrand()<<"\t0"<<endl;
}



void Application::emitIndels(const CigarAlignment &revAlign)
{
  const int L=revAlign.length();
  bool inInsertion=false;
  int beginInsertion=0;
  int prevRefPos=-1, prevAltPos=-1;
  for(int i=0 ; i<L ; ++i) {
    const int refPos=revAlign[i];
    if(refPos==CIGAR_UNDEFINED) {
      if(!inInsertion) { beginInsertion=i; inInsertion=true; }
    }
    else {
      if(inInsertion) {
	cout<<altSubstrate<<"\tindel\tinsertion\t"<<beginInsertion+1
	    <<"\t"<<i<<"\t.\t+\t0"<<endl;
	inInsertion=false;
      }
      if(prevRefPos>=0 && refPos-prevRefPos>1)
	cout<<altSubstrate<<"\tindel\tdeletion\t"<<prevAltPos+1
	    <<"\t"<<i<<"\t.\t+\t0"<<endl;
      prevRefPos=refPos;
      prevAltPos=i;
    }
  }
}



int Application::mapConsensus(int refPos,const CigarAlignment &alignment)
{
  while(refPos>0 && alignment[refPos]==CIGAR_UNDEFINED) --refPos;
  return alignment[refPos];
}


