/****************************************************************
 extract-splice-features.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Interval.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/FastaReader.H"
#include "IsochoreTable.H"
#include "SpliceFeatureExtractor.H"
using namespace std;
using namespace BOOM;

Alphabet alphabet=DnaAlphabet::global();

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  float getGcContent(const String &);
};


int main(int argc,char *argv[])
{
  try {
    Application app;
    return app.main(argc,argv);
  }
  catch(const char *p) { cerr << p << endl; }
  catch(const string &msg) { cerr << msg.c_str() << endl; }
  catch(const exception &e)
    {cerr << "STL exception caught in main:\n" << e.what() << endl;}
  catch(...) { cerr << "Unknown exception caught in main" << endl; }
  return -1;
}



Application::Application()
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"d:");
  if(cmd.numArgs()!=6)
    throw String("\n\
extract-splice-features <*.iso> <*.fasta> <GT|AG> <begin:end> <motifs.txt> <motif-distance-param>\n\
  <*.iso> = GeneZilla parameter file\n\
  <begin:end> = interval to scan; 0-based coords, half open: [begin,end)\n\
  <motif-distance-param> must be strictly between 0 and 1\n\
  Motifs file should contain two columns: the motif and its score\n\
  -d N : ignore regulatory elements further than N bases from splice site\n\
  Output: site location, site score, regulatory score\n\
          (one line per putative site)\n\
");
  const String isoFile=cmd.arg(0);
  const String fastaFile=cmd.arg(1);
  const String GTorAG=cmd.arg(2);
  const String intervalStr=cmd.arg(3);
  const String motifFile=cmd.arg(4);
  const float distanceParm=cmd.arg(5).asFloat();
  Vector<String> fields; intervalStr.getFields(fields,":");
  if(fields.size()!=2) throw intervalStr+" : invalid interval specification";
  const Interval interval(fields[0].asInt(),fields[1].asInt());
  const int maxDistance=cmd.option('d') ? cmd.optParm('d').asInt() : 1000000;

  // Load sequence
  String defline, seqStr;
  FastaReader::load(fastaFile,defline,seqStr);
  Sequence seq(seqStr,alphabet);
  float gcContent=getGcContent(seqStr);

  // Load models
  GarbageCollector gc;
  IsochoreTable isochores(gc);
  isochores.load(isoFile);
  Isochore *isochore=isochores.getIsochore(gcContent);
  SignalSensor *sensor;
  SignalType signalType=stringToSignalType(GTorAG);
  sensor=isochore->signalTypeToSensor[signalType];

  // Load regulatory motifs
  RegulatoryMotifs motifs(motifFile);

  // Invoke the feature extractor
  SpliceFeatureExtractor extractor(*sensor,motifs,distanceParm,maxDistance);
  const int L=seqStr.length();
  const int consensusOffset=sensor->getConsensusOffset();
  const int contextWindowLen=sensor->getContextWindowLength();
  const int begin=consensusOffset, end=L-contextWindowLen+consensusOffset;
  for(int pos=begin ; pos<=end ; ++pos)
    if(seqStr.substring(pos,2)==GTorAG) {
      float intrinsicScore, regulatoryScore;
      if(extractor.extract(seq,seqStr,pos,intrinsicScore,regulatoryScore))
	cout<<GTorAG<<"\t"<<pos<<"\t"<<intrinsicScore<<"\t"<<regulatoryScore<<endl;
    }

  return 0;
}



float Application::getGcContent(const String &s)
{
  const int A=s.count('A')+s.count('a');
  const int C=s.count('C')+s.count('c');
  const int G=s.count('G')+s.count('g');
  const int T=s.count('T')+s.count('t');
  const int total=A+C+G+T;
  const float gc=(G+C)/float(total);
  return gc;
}




