/****************************************************************
 vcf-to-gcf.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/File.H"
#include "BOOM/Pipe.H"
#include "BOOM/Vector.H"
#include "BOOM/Map.H"
#include "BOOM/Regex.H"
#include "BOOM/GCF.H"
#include "BOOM/Time.H"
using namespace std;
using namespace BOOM;

struct Region {
  int begin, end;
  Region(int b,int e) : begin(b), end(e) {}
  bool contains(int p) { return begin<=p && p<end; }
};

struct Variant {
  String chr;
  int pos;
  String ref, alt, id;
  Variant(const String &chr,int pos,const String &ref,
	  const String &alt,const String &id)
    : chr(chr), pos(pos), ref(ref), alt(alt), id(id) {}
  String asString() { return id+":"+chr+":"+pos+":"+ref+":"+alt; }
};

struct Individual {
  String id;
  Vector<String> genotypes;
};

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  Time timer;
  Map<String,Vector<Region> > regions;
  Vector<Individual> individuals;
  Vector<Variant> variants;
  Regex gzipRegex; // *.gz
  Regex dnaRegex;
  bool wantFilter;
  bool prependChr;
  bool SNPsOnly;
  bool variableOnly;
  bool quiet;
  void loadRegions(const String &filename);
  void convert(File &infile,File &outfile);
  void convertSM(File &infile,File &outfile,const String &tempfile);
  void preprocess(File &infile);
  void parseChromLine(const Vector<String> &);
  bool parseVariant(const Vector<String> &fields,String &chr,int &pos,
		    String &ref,String &alt,String &id);
  bool parseVariant(const Vector<String> &);
  void parseVariantSM(const Vector<String> &fields,File &temp,
		      int &variantNum,int entrySize);
  void parseVariantAndGenotypes(const Vector<String> &);
  bool keep(const String &chr,int pos);
  void output(File &outfile);
  void outputSM(File &temp,File &out,int entrySize);
  bool variableSite(const Vector<String> &fields);
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
  catch(...)
    {cerr << "Unknown exception caught in main" << endl;}
  return -1;
}



Application::Application()
  : gzipRegex(".*\\.gz"), dnaRegex("^[ACGTacgt]+$")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"cf:qsvm:");
  if(cmd.numArgs()!=2)
    throw String("\nvcf-to-gcf [options] <in.vcf> <out.gcf>\n\
   both input and output files can be zipped (use .gz as suffix)\n\
   -f regions.bed : keep only variants in these regions\n\
        (coordinates are 0-based, end is not inclusive)\n\
   -c : prepend \"chr\" before chromosome names\n\
   -s : SNPs only\n\
   -v : variable sites only\n\
   -q : quiet (no warnings)\n\
   -m <tempfile> : small memory footprint (may be slow)\n\
");
  const String infile=cmd.arg(0);
  const String outfile=cmd.arg(1);
  variableOnly=cmd.option('v');
  wantFilter=cmd.option('f');
  prependChr=cmd.option('c');
  SNPsOnly=cmd.option('s');
  quiet=cmd.option('q');
  const bool smallmem=cmd.option('m');

  // Load regions to filter by
  if(wantFilter) loadRegions(cmd.optParm('f'));

  // Open files
  File *vcf=gzipRegex.match(infile) ? 
    new Pipe(String("cat ")+infile+" | gunzip","r") : 
    new File(infile);
  File *gcf=gzipRegex.match(outfile) ?
    new Pipe(String("bgzip > ")+outfile,"w") :
    new File(outfile,"w");

  // Perform conversion
  timer.startCounting();
  if(smallmem) {
    cerr<<"preprocessing..."<<endl;
    preprocess(*vcf);
    delete vcf;
    vcf=gzipRegex.match(infile) ?  
      new Pipe(String("cat ")+infile+" | gunzip","r") :
      new File(infile); 
    convertSM(*vcf,*gcf,cmd.optParm('m'));
  }
  else convert(*vcf,*gcf);
  vcf->close(); gcf->close();
  cerr<<"done."<<endl;
  cerr<<timer.elapsedTime()<<endl;

  return 0;
}



void Application::loadRegions(const String &filename)
{
  File file(filename);
  while(!file.eof()) {
    String line=file.getline();
    line.trimWhitespace();
    Vector<String> &fields=*line.getFields();
    if(fields.size()>0) {
      if(fields.size()<3) throw filename+" - can't parse bed file: "+line;
      const String &chr=fields[0];
      const int begin=fields[1].asInt();
      const int end=fields[2].asInt();
      if(!regions.isDefined(chr)) regions[chr]=Vector<Region>();
      regions[chr].push_back(Region(begin,end));
    }
    delete &fields;
  }
}



void Application::convert(File &infile,File &outfile)
{
  while(!infile.eof()) {
    String line=infile.getline();
    line.trimWhitespace();
    Vector<String> &fields=*line.getFields();
    if(fields.size()>0) {
      if(fields[0]=="#CHROM") parseChromLine(fields);
      else if(fields[0][0]!='#') parseVariantAndGenotypes(fields);
    }
    delete &fields;
  }
  output(outfile);
}



void Application::parseChromLine(const Vector<String> &fields)
{
  // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  200848168@1097100704
  const int numFields=fields.size();
  if(numFields<10 || fields[8]!="FORMAT") 
    throw "Error parsing #CHROM line in VCF file";
  int numIndiv=numFields-9;
  individuals.resize(numIndiv);
  for(int i=0 ; i<numIndiv ; ++i)
    individuals[i].id=fields[i+9];
}



void Application::parseVariantAndGenotypes(const Vector<String> &fields)
{
  if(!parseVariant(fields)) return;
  const int numIndiv=fields.size()-9;
  for(int i=0 ; i<numIndiv ; ++i) {
    const String &genotype=fields[i+9];
    Individual &indiv=individuals[i];
    indiv.genotypes.push_back(genotype);
  }
}



bool Application::keep(const String &chr,int pos)
{
  if(!regions.isDefined(chr)) return false;
  const Vector<Region> &regs=regions[chr];
  for(Vector<Region>::const_iterator cur=regs.begin(), end=regs.end() ;
      cur!=end ; ++cur) {
    if((*cur).contains(pos)) return true;
  }
  return false;
}



void Application::output(File &out)
{
  // Output list of variants
  const int numVariants=variants.size();
  for(int i=0 ; i<numVariants ; ++i) {
    const Variant &v=variants[i];
    out.print(v.asString());
    out.print(i+1<numVariants ? "\t" : "\n");
  }
  if(numVariants==0) out.print("\n");

  // Output individual genotypes
  for(Vector<Individual>::const_iterator cur=individuals.begin(), 
	end=individuals.end() ; cur!=end ; ++cur) {
    const Individual &indiv=*cur;
    if(numVariants==0) { out.print(indiv.id+"\n"); continue; }
    out.print(indiv.id+"\t");
    const int numVariants=indiv.genotypes.size();
    for(int i=0 ; i<numVariants ; ++i) {
      out.print(indiv.genotypes[i]);
      out.print(i+1<numVariants ? "\t" : "\n");
    }
  }
}



void Application::preprocess(File &infile)
{
  while(!infile.eof()) {
    String line=infile.getline();
    line.trimWhitespace();
    Vector<String> &fields=*line.getFields();
    if(fields.size()>0) {
      if(fields[0]=="#CHROM") parseChromLine(fields);
      else if(fields[0][0]!='#') parseVariant(fields);
    }
    delete &fields;
  }
}



bool Application::variableSite(const Vector<String> &fields)
{
  bool seenZero=false, seenOne=false;
  for(int i=9 ; i<fields.size() ; ++i) {
    const String &field=fields[i];
    if(field.length()==3) {
      if(field[0]=='0' || field[3]=='0') seenZero=true;
      if(field[0]=='1' || field[3]=='1') seenOne=true;
    }
    else if(field.length()==1) {
      if(field[0]=='0') seenZero=true;
      if(field[0]=='1') seenOne=true;
    }
    else throw String("Bad field in VCF: ")+field;
  }
  return seenZero && seenOne;
}



bool Application::parseVariant(const Vector<String> &fields,
			       String &chr,int &pos,String &ref,
			       String &alt,String &id)
{
  if(variableOnly && !variableSite(fields)) return false;
  if(fields.size()<10 || fields[6]!="PASS" || fields[8]!="GT") return false;
  chr=fields[0];
  if(prependChr) chr=String("chr")+chr;
  pos=fields[1].asInt()-1; // VCF files are 1-based
  if(wantFilter && !keep(chr,pos)) return false;
  id=fields[2];
  if(id==".") id=chr+"@"+pos;
  ref=fields[3];
  alt=fields[4];
  if(SNPsOnly && (ref.size()!=1 || alt.size()!=1)) return false;
  if(!dnaRegex.match(alt) || !dnaRegex.match(ref)) 
    return false; // nonstandard or structural variant
  return true;
}



bool Application::parseVariant(const Vector<String> &fields)
{
  String chr, ref, alt, id;
  int pos;
  if(!parseVariant(fields,chr,pos,ref,alt,id)) return false;
  variants.push_back(Variant(chr,pos,ref,alt,id));
  return true;
}



void Application::parseVariantSM(const Vector<String> &fields,File &temp,
				 int &variantNum,int entrySize)
{
  // Parse the line
  String chr, ref, alt, id;
  int pos;
  if(!parseVariant(fields,chr,pos,ref,alt,id)) return;

  const int numVariants=variants.size();
  if(variantNum>=numVariants) {
    cerr<<variantNum<<" >= "<<numVariants<<endl;
    INTERNAL_ERROR;
  }

  // Parse genotypes & store in binary temp file
  const int rowSize=numVariants*entrySize;
  char *buffer=new char[entrySize];
  const int numIndiv=individuals.size();
  for(int i=0 ; i<numIndiv ; ++i) {
    const String &genotype=fields[i+9];
    for(int j=0 ; j<entrySize ; ++j) buffer[j]=' ';
    if(genotype.length()>entrySize) {
      cerr<<"Genotype length "<<genotype.length()<<" > "<<entrySize<<endl;
      INTERNAL_ERROR;
    }
    for(int j=0 ; j<genotype.length() ; ++j) buffer[j]=genotype[j];
    temp.seek(i*rowSize+variantNum*entrySize);
    temp.write(entrySize,reinterpret_cast<void*>(buffer));
  }
  delete [] buffer;
  ++variantNum;
}



void Application::convertSM(File &infile,File &outfile,const String &tempfile)
{
  cerr<<timer.elapsedTime()<<endl;
  cerr<<"writing binary file..."<<endl;

  // First, allocate binary file to store genotypes
  File temp(tempfile,"w");
  const int numIndiv=individuals.size(), numVariants=variants.size();
  const int entrySize=3;
  const int rowSize=numVariants*entrySize;
  const int totalSize=numIndiv*rowSize;
  const int lastEntry=totalSize-entrySize;
  temp.seek(lastEntry);
  char *buffer=new char[entrySize];
  temp.write(entrySize,reinterpret_cast<void*>(buffer));
  delete [] buffer;

  // Reprocess the input file, store genotypes in the binary file
  int variantNum=0;
  while(!infile.eof()) {
    String line=infile.getline();
    line.trimWhitespace();
    Vector<String> &fields=*line.getFields();
    if(fields.size()>0 && fields[0][0]!='#') 
      parseVariantSM(fields,temp,variantNum,entrySize);
    delete &fields;
  }
  
  // Convert the binary temp file into output GCF file
  temp.close();
  cerr<<timer.elapsedTime()<<endl;
  cerr<<"converting binary file to text..."<<endl;
  File tempRead(tempfile,"r");
  outputSM(tempRead,outfile,entrySize);
}



void Application::outputSM(File &temp,File &out,int entrySize)
{
  // Output list of variants
  const int numVariants=variants.size();
  for(int i=0 ; i<numVariants ; ++i) {
    const Variant &v=variants[i];
    out.print(v.asString());
    out.print(i+1<numVariants ? "\t" : "\n");
  }
  if(numVariants==0) out.print("\n");

  // Output individual genotypes
  const int numIndiv=individuals.size();
  temp.seek(0);
  char *buffer=new char[entrySize+1];
  buffer[entrySize]='\0';
  for(int i=0 ; i<numIndiv ; ++i) {
    const Individual &indiv=individuals[i];
    if(numVariants==0) { out.print(indiv.id+"\n"); continue; }
    out.print(indiv.id+"\t");
    for(int i=0 ; i<numVariants ; ++i) {
      temp.read(entrySize,reinterpret_cast<void*>(buffer));
      String genotype(buffer);
      genotype.trimWhitespace();
      out.print(genotype);
      out.print(i+1<numVariants ? "\t" : "\n");
    }
  }
  delete buffer;
}


