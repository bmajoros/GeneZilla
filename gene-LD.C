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
#include "BOOM/Map.H"
using namespace std;
using namespace BOOM;

struct Individual {
  int hap[2][2]; // [which chrom][which gene]
};

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
  void recode(GCF &,int whichGene);
protected:
  Map<String,int> haploMap[2]; // for two genes
  Vector<String> haplotypes[2]; // for two genes
  Array1D<Individual> individuals;
  String genotypeToString(const Array1D<int> &genotype);
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
  individuals.resize(numIndiv);

  // Re-code the haplotypes as integers (separately for each gene)
  recode(gcf1,0); recode(gcf2,1);

  // Compute mutual information


  return 0;
}



void Application::recode(GCF &gcf,int whichGene)
{
  const int numIndiv=gcf.numIndividuals();
  for(int i=0 ; i<numIndiv ; ++i) {
    const GcfIndividual &indiv=gcf.getIthIndividual(i);
    recode(indiv.chrom[0],whichGene);
    recode(indiv.chrom[1],whichGene);
  }
}



void Application::recode(const Array1D<int> &chrom,int whichGene)
{
  String s=genotypeToString(indiv.chrom[0]);
  int hapID;
  if(!haploMap[whichGene].isDefined(s)) {
    hapID=haploMap[whichGene][s]=haplotypes[whichGene].size();
    haplotypes[whichGene].push_back(s);
  }
  else hapID=haploMap[whichGene][s];
  
}



String Application::genotypeToString(const Array1D<int> &gt)
{
  const int n=gt.size();
  char str[n+1]; str[n]='\0';
  for(int i=0 ; i<n ; ++i) str[i]='0'+gt[i];
  return str;
}


