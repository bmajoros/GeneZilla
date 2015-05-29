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
using namespace std;

#ifdef EXPLICIT_GRAPHS
//#error Please edit the file genezilla.H, comment out the definition of EXPLICIT_GRAPHS, issue a "make clean", and recompile this project
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
};



/***************************************************************************
                            class Application
 ***************************************************************************/
class Application {
public:
  void AppMain(int argc,char *argv[]);
protected:
  void parseVariantLine(const String &line,VariantRegion &);
  void applySensor(SignalSensor &,const VariantRegion &,Sequence &refSeq,
		   Sequence &altSeq,const String &refSeqStr,const String 
		   &altSeqStr,CigarAlignment &);
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
  GeneZilla genezilla(PROGRAM_NAME,VERSION,new EdgeFactory,-1);
  BOOM::String defline, refSeqStr, altSeqStr;
  FastaReader::load(refFasta,defline,refSeqStr);
  FastaReader::load(altFasta,defline,altSeqStr);
  Sequence refSeq(refSeqStr,alphabet);
  Sequence altSeq(altSeqStr,alphabet);
  const float gc=genezilla.getGCcontent(seqString);
  Isochore *isochore=genezilla.getIsochore(gc);
  CigarString cigar(cigarFile);
  CigarAlignment &alignment=*cigar.getAlignment();

  // Iterate through variants
  File variantFile(variantFilename);
  VariantRegion region;
  while(!variantFile.eof()) {
    String line=variantFile.getline();
    parseVariantLine(line,region);

    // Iterate through signal sensors
    for(Vector<SignalSensor*>::iterator cur=signalSensors.begin(), end=
	  signalSensors.end() ; cur!=end ; ++cur) {
      SignalSensor *sensor=*cur;
      applySensor(*sensor,region,refSeq,altSeq,refSeqStr,altSeqStr,
		  alignment);


    }
  }
}



void Application::parseVariantLine(const String &line,VariantRegion &region)
{
  line.trimWhitespace();
  Vector<String> &fields=*line.getFields();
  if(fields.size()==6) {
    region.refBegin=fields[1].asInt();
    region.refEnd=fields[2].asInt();
    region.altBegin=fields[4].asInt();
    region.altEnd=fields[5].asInt();
  }
  delete &fields;
}



void Application::applySensor(SignalSensor &sensor,const VariantRegion &region,
			      Sequence &refSeq,Sequence &altSeq,const String
			      &refSeqStr,const String &altSeqStr,
			      CigarAlignment &alignment)
{
  const float threshold=sensor.getCutoff();
  const int windowLen=sensor.getContextWindowLength();
  const int refLen=refSeq.getLength(), altLen=altSeq.getLength();

  // Find reference signals missing in the alternate
  int begin=region.refBegin-windowLen+1;
  if(begin<0) begin=0;
  int end=regin.refEnd-1;
  if(end+windowLen>refLen) end=refLen-windowLen;
  for(int r=begin ; r<=end ; ++r) {
    SignalPtr signal=sensor.detect(refSeq,refSeqStr,r);
    if(signal) {
      if(containsIndel(signal,alignment)) {
	// emit: this signal is missing in the alt
	
      }
      else {
	// run on corresponding alt interval, see if score < threshold

      }
    }
  }

  // Find alternate signals missing in the reference



  // Find indels


}



