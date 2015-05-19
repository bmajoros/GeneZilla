/****************************************************************
 split-gff-by-gene.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
using namespace std;
using namespace BOOM;

class Application {
  bool wantLongest;
public:
  Application();
  int main(int argc,char *argv[]);
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
  CommandLine cmd(argc,argv,"l");
  if(cmd.numArgs()!=2)
    throw String("split-gff-by-gene [options] <in.gff> <out-dir>");
  const String inGff=cmd.arg(0);
  const String outDir=cmd.arg(1);
  wantLongest=cmd.option('l');

  // Load GFF file
  Vector<GffGene> &genes=*GffReader::loadGenes(inGff);


  return 0;
}

