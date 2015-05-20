/****************************************************************
 project-annotation.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/CigarString.H"
#include "BOOM/GffReader.H"
#include "BOOM/FastaReader.H"
#include "Labeling.H"
using namespace std;
using namespace BOOM;


class Application
{
public:
  Application();
  int main(int argc,char *argv[]);
private:
  String loadSeq(const String &filename);
  GeneModelLabel getExonLabel(int phase);
  void computeLabeling(TranscriptList *,Labeling &);
  void mapLabeling(Labeling &from,Labeling &to,const String &cigar);
  void mapTranscripts(TranscriptList &,const String &cig,const String &outfile);
};



int main(int argc,char *argv[])
  {
    try
      {
	Application app;
	return app.main(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
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
  if(cmd.numArgs()!=6)
    throw String("project-annotation <ref.gff> <ref.fasta> <alt.fasta> <ref-alt.cigar> <out.vector> <out.gff>");
  const String refGffFile=cmd.arg(0);
  const String refFasta=cmd.arg(1);
  const String altFasta=cmd.arg(2);
  const String cigarFile=cmd.arg(3);
  const String outfile=cmd.arg(4);
  const String outGff=cmd.arg(5);
  
  // Read some data from files
  String refSeq=loadSeq(refFasta), altSeq=loadSeq(altFasta);
  int refSeqLen=refSeq.length(), altSeqLen=altSeq.length();
  GffReader gffReader(refGffFile);
  TranscriptList *transcripts=gffReader.loadTranscripts();

  // Compute the reference labeling
  Labeling refLab(refSeqLen);
  computeLabeling(transcripts,refLab);

  // Load CIGAR string
  String CIGAR;
  File cigar(cigarFile);
  CIGAR=cigar.getline();
  cigar.close();

  // Project the reference labeling over to the alternate sequence
  Labeling altLab(altSeqLen);
  mapLabeling(refLab,altLab,CIGAR);

  // Project the reference GFF over to an alternate GFF
  mapTranscripts(*transcripts,CIGAR,outGff);

  // Generate output
  ofstream os(outfile.c_str());
  os<<altLab;

  return 0;
}



void Application::mapTranscripts(TranscriptList &transcripts,const String &cig,
				 const String &outfile)
{
  CigarString cigar(cig);
  CigarAlignment &align=*cigar.getAlignment();
  GffTranscript *transcript=transcripts[0];
  char strand=transcript->getStrand();
  String source=transcript->getSource();
  String substrate=transcript->getSubstrate();
  int numExons=transcript->getNumExons();
  for(int i=0 ; i<numExons ; ++i) {
    GffExon &exon=transcript->getIthExon(i);
    int begin=exon.getBegin(), end=exon.getEnd();
    begin=align[begin]; end=align[end-1]+1;//###
    exon.setBegin(begin); exon.setEnd(end);
  }
  delete &align;
  ofstream os(outfile.c_str());
  transcript->toGff(os);
}



void Application::mapLabeling(Labeling &from,Labeling &to,const String &cig)
{
  CigarString cigar(cig);
  CigarAlignment &align=*cigar.getAlignment();
  to.asArray().setAllTo(LABEL_NONE);
  int L=align.length();
  for(int i=0 ; i<L ; ++i) {
    int j=align[i];
    //cout<<"j="<<j<<" toLen="<<to.length()<<endl;
    if(j!=CIGAR_UNDEFINED) to[j]=from[i];
  }
  delete &align;
}



void Application::computeLabeling(TranscriptList *transcripts,
				  Labeling &refLab)
{
  int numTrans=transcripts->size();
  if(numTrans!=1) throw "Number of transcripts provided must be exactly 1";
  GffTranscript *transcript=(*transcripts)[0];
  int begin=transcript->getBegin(), end=transcript->getEnd();
  char strand=transcript->getStrand();
  if(strand!='+') throw "only forward-strand features are currently supported";
  int numExons=transcript->getNumExons();
  refLab.asArray().setAllTo(LABEL_INTERGENIC);
  for(int i=begin ; i<end ; ++i) refLab[i]=LABEL_INTRON;
  int phase=0;
  for(int i=0 ; i<numExons ; ++i) {
    GffExon &exon=transcript->getIthExon(i);
    begin=exon.getBegin(); end=exon.getEnd();
    for(int j=begin ; j<end ; ++j) {
      refLab[j]=getExonLabel(phase);
      phase=(phase+1)%3;
    }
  }
}



GeneModelLabel Application::getExonLabel(int phase)
{
  switch(phase)
    {
    case 0: return LABEL_EXON_0;
    case 1: return LABEL_EXON_1;
    case 2: return LABEL_EXON_2;
    default: throw String("bad phase in getExonLabel(): ")+phase;
    }
}



String Application::loadSeq(const String &filename)
{
  FastaReader reader(filename);
  String def, seq;
  if(!reader.nextSequence(def,seq)) throw filename+" : cannot read file";
  return seq;
}

