/****************************************************************
 gene-LD.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/GCF.H"
using namespace std;
using namespace BOOM;

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
  void recode(GCF &,int whichGene);
projected:
  Map<String,int> haploMap[2]; // for two genes
  Vector<String> haplotypes[2]; // for two genes
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
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String("gene-LD <gene1.gcf> <gene2.gcf>");
  const String file1=cmd.arg(0);
  const String file2=cmd.arg(1);

  // Load files
  GCF gcf1(file1), gcf2(file2);
  const int numIndiv=gcf1.numIndividuals();
  if(gcf2.numIndividuals()!=numIndiv) throw "files have different numbers of individuals";

  // Re-code the haplotypes as integers (separately for each gene)
  recode(gcf1,0); recode(gcf2,1);

  // Compute mutual information


  return 0;
}



void Application::recode(GCF &gcf,int whichGene)
{
  const int numIndiv=gcf.numIndividuals();
  
}


