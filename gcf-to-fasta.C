/****************************************************************
 gcf-to-fasta.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/FastaWriter.H"
#include "BOOM/File.H"
#include "BOOM/Vector.H"
using namespace std;
using namespace BOOM;

struct Region {
  String chr;
  int begin, end;
  String seq;
  Region(const String &chr,int begin,int end,const String &seq)
    : chr(chr), begin(begin), end(end), seq(seq) {}
};

struct Genotype {
  Array1D<int> alleles;
  Genotype(int ploidy) : alleles(ploidy) {}
};

struct Variant {
  String id, chr, ref, alt;
  int pos;
  Variant(const String &id,const String &chr,int pos,
	  const String &ref,const String &alt)
    : id(id), chr(chr), pos(pos), ref(ref), alt(alt) {}
};

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  String twoBitToFa;
  Vector<Variant> variants;
  FastaWriter writer;
  Vector<Region> regions;
  void convert(File &gcf,ostream &,const String genomeFile);
  void parseHeader(const String &line);
  void loadRegions(const String &regionsFilename,const String &genomeFilename,
		   const String &tempFile);
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
  : twoBitToFa("twoBitToFa")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"t:");
  if(cmd.numArgs()!=4)
    throw String("\ngcf-to-fasta [options] <in.gcf> <genome.2bit> <regions.bed> <out.fasta>\n\
     -t path : path to twoBitToFa\n\
");
  const String &gcfFilename=cmd.arg(0);
  const String &genomeFilename=cmd.arg(1);
  const String &regionsFilename=cmd.arg(2);
  const String &fastaFilename=cmd.arg(3);
  if(cmd.option('t')) twoBitToFa=cmd.optParm('t');

  // Load regions
  loadRegions(regionsFilename,genomeFilename,fastaFilename);

  // Process GCF file
  File gcf(gcfFilename);
  ofstream os(fastaFilename.c_str());
  convert(gcf,os,genomeFilename);

  return 0;
}



void Application::convert(File &gcf,ostream &os,const String genomeFile)
{
  // Parse header
  String line=gcf.getline();
  line.trimWhitespace();
  parseHeader(line);
  const int numVariants=variants.size();

  // Process each individual
  while(!gcf.eof()) {
    line=gcf.getline();
    line.trimWhitespace();
    if(line.isEmpty()) continue;
    Vector<String> &fields=*line.getFields();
    if(fields.size()!=numVariants+1) throw "Wrong number of fields in GCF line";
    String id=fields.front();
    fields.erase(fields.begin());
    for(Vector<String>::const_iterator cur=fields.begin(), end=fields.end() ;
	cur!=end ; ++cur) {
      const String &field=*cur;
      if(field.length()!=3) throw String("Cannot parse genotype: ")+field;
      int gt[2]; gt[0]=field[0]-'0'; gt[1]=field[2]-'0';
    }
    delete &fields;
    emit(
  }
}



void Application::parseHeader(const String &line)
{
  Vector<String> &fields=*line.getFields();
  for(Vector<String>::const_iterator cur=fields.begin(), end=fields.end()
	; cur!=end ; ++cur) {
    const String &field=*cur;
    Vector<String> &fields=*field.getFields(":");
    if(fields.size()!=5) throw String("cannot parse GCF header: ")+field;
    const String &id=fields[0];
    const String &chr=fields[1];
    const int pos=fields[2];
    const String &ref=fields[3];
    const String &alt=fields[4];
    variants.push_back(Variant(id,chr,pos,ref,alt));
    delete &fields;
  }
  delete &fields;
}



void Application::loadRegions(const String &regionsFilename,const String &
			      genomeFilename,const String &tempFile)
{
  File reg(regionsFilename);
  while(!reg.eof()) {
    String line=reg.getline();
    line.trimWhitespace();
    if(line.isEmpty()) continue;
    Vector<String> &fields=*line.getFields();
    if(fields.size()<3) throw String("error in regions file: ")+line;
    const String chr=fields[0];
    const int begin=fields[1].asInt(), end=fields[2].asInt();
    delete &fields;
    
    // Invoke twoBitToFa to extract sequence from chrom file
    String cmd=twoBitToFa+" -seq="+chr+" -start="+begin+" -end="+end+
      +" "+genomeFilename+" "+tempFile;
    system(cmd.c_str());
    String def, seq;
    FastaReader::load(tempFile,def,seq);
    regions.push_back(Region(chr,begin,end,seq));
  }
  unlink(tempFile.c_str());
}



