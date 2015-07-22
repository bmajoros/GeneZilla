/****************************************************************
 remove-unpaired-reads.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Regex.H"
#include "BOOM/File.H"
using namespace std;
using namespace BOOM;

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  Regex idRegex;
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
  : idRegex("(.*)(\\/\d)")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=6)
    throw String("remove-unpaired-reads <in1.fastq.gz> <in2.fastq.gz> <out1.fastq.gz> <out2.fastq.gz> <unpaired1.fastq.gz> <unpaired2.fastq.gz>");
  const String infile1=cmd.arg(0);
  const String infile2=cmd.arg(1);
  const String outfile1=cmd.arg(2);
  const String outfile2=cmd.arg(3);
  const String unpaired1=cmd.arg(4);
  const String unpaired2=cmd.arg(5);

  GunzipPipe in1(infile1), in2(infile2);
  GzipPipe out1(outfile1), out2(outfile2);
  GzipPipe singletons1(unpaired1), singletons2(unpaired2);

  return 0;
}

