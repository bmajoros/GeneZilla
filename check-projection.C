/****************************************************************
 check-projection.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/GffReader.H"
#include "BOOM/ProteinTrans.H"
#include "BOOM/CodonIterator.H"
#include "Labeling.H"
#include "ProjectionChecker.H"
using namespace std;
using namespace BOOM;



class Application {
public:
  int main(int argc,char *argv[]);
protected:
  GffTranscript *loadGff(const String &filename);
  void parseNoncanonicals(const String &,Set<String> &);
};




int main(int argc,char *argv[]) {
  try {
    Application app;
    return app.main(argc,argv); }
  catch(const char *p) { cerr << p << endl; }
  catch(const string &msg) { cerr << msg.c_str() << endl; }
  catch(const exception &e) 
    {cerr << "STL exception caught in main:\n" << e.what() << endl;}
  catch(...)
    {cerr << "Unknown exception caught in main" << endl;}
  return -1;
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"d:a:");
  if(cmd.numArgs()!=5)
    throw String("\n\
check-projection <ref.fasta> <ref.gff> <alt.fasta> <projected.gff> <labels.txt>\n\
     -d donors = allow noncanonical donors; comma-separated list (capital letters)\n\
     -a acceptors = allow noncanonical acceptors; comma-separated list (capital letters)\n\
");
  const String refFasta=cmd.arg(0);
  const String refGff=cmd.arg(1);
  const String altFasta=cmd.arg(2);
  const String altGff=cmd.arg(3);
  const String labelFile=cmd.arg(4);
  Set<String> nonCanonicalGTs, nonCanonicalAGs;
  if(cmd.option('d')) parseNoncanonicals(cmd.optParm('d'),nonCanonicalGTs);
  if(cmd.option('a')) parseNoncanonicals(cmd.optParm('a'),nonCanonicalAGs);

  // Load input files
  String def, refSubstrate, altSubstrate;
  FastaReader::load(refFasta,def,refSubstrate);
  FastaReader::load(altFasta,def,altSubstrate);
  GffTranscript *refTrans=loadGff(refGff), *altTrans=loadGff(altGff);
  Labeling labeling(labelFile);

  // Create the checker object
  ProjectionChecker checker(*refTrans,*altTrans,refSubstrate,
			    altSubstrate,labeling,
			    nonCanonicalGTs,nonCanonicalAGs);

  // Check splice sites
  bool splicingIsOK=checker.checkSpliceSites(false);

  // Translate to proteins
  String refProtein, altProtein;
  checker.translate(*refTrans,*altTrans,refProtein,altProtein);
  
  // Check for start codon
  if(!checker.hasStartCodon(altProtein)) cout<<"No start codon"<<endl;

  // Everything else depends on splice sites being intact
  if(!splicingIsOK) return 0;

  // Check for frameshifts
  if(refProtein!=altProtein) cout<<"proteins differ"<<endl;
  checker.checkFrameshifts(labeling,*altTrans,altSubstrate,false);

  // Check for stop codons
  int PTCpos;
  if(checker.hasPTC(altProtein,PTCpos)) {
    cout<<"premature stop at AA position "<<PTCpos<<" in alt protein"<<endl;
    if(!checker.detectNMD(*altTrans,altSubstrate,false))
      cout<<"truncation predicted"<<endl;
  }
  else if(!checker.hasStopCodon(altProtein))
    cout<<"missing stop codon"<<endl; 
  
  return 0;
}



GffTranscript *Application::loadGff(const String &filename)
{
  GffReader reader(filename);
  Vector<GffTranscript*> *transcripts=reader.loadTranscripts();
  const int n=transcripts->size();
  if(n<1) throw filename+" contains no transcripts";
  GffTranscript *transcript=(*transcripts)[0];
  for(int i=1 ; i<n ; ++i) delete (*transcripts)[i];
  delete transcripts;
  return transcript;
}



void Application::parseNoncanonicals(const String &str,Set<String> &into)
{
  Vector<String> &fields=*str.getFields(",");
  for(Vector<String>::const_iterator cur=fields.begin(), end=fields.end() ;
      cur!=end ; ++cur)
    into.insert(*cur);
  delete &fields;
}



