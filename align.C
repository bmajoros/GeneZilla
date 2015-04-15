#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "AlignmentSubstMatrix.H"
#include "BOOM/AminoAlphabet.H"
#include "BOOM/DnaAlphabet.H"
#include "Needleman.H"


class Application {
  bool wantUngapped;
public:
  Application();
  int main(int argc,char *argv[]);
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
    CommandLine cmd(argc,argv,"qn:w:uUlc:");
    if(cmd.numArgs()!=6)
      throw string("align <AlignmentSubstMatrix> <+GapOpenPenalty> <+GapExtendPenalty> \n               <*.fasta> <*.fasta> DNA|PROTEIN [-q] [-w #]\n\n\
example: align blosum62 5 2 1.fasta 2.fasta DNA\n\
-c <file> = emit CIGAR string into file\n\
-q = quiet (no alignment output -- just match intervals, #matches, and score)\n\
-u = ungapped alignment\n\
-U = ungapped, report only the closed match interval and #matches\n\
-l = local alignment (no penalty for gaps at ends)\n\
-w # = format output to given width (default is 60)\n");
    String matrixFile=cmd.arg(0);
    double gapOpen=-fabs(cmd.arg(1).asDouble());
    double gapExtend=-fabs(cmd.arg(2).asDouble());
    String file1=cmd.arg(3);
    String file2=cmd.arg(4);
    String type=cmd.arg(5);
    wantUngapped=cmd.option('u') || cmd.option('U');
    bool terseUngapped=cmd.option('U');
    bool local=cmd.option('l');
    if(local && wantUngapped) throw "local ungapped not currently supported";
    if(type!="DNA" && type!="PROTEIN") throw "specify DNA or PROTEIN";
    bool quiet=cmd.option('q');
    if(cmd.option('w')) AlignmentPath::MAX_WIDTH=cmd.optParm('w').asInt();
    ofstream cigarFile;
    if(cmd.option('c')) cigarFile.open(cmd.optParm('c').c_str());

    Alphabet &alphabet=
      (type=="DNA") ? 
      static_cast<Alphabet&>(DnaAlphabet::global()) : 
      static_cast<Alphabet&>(AminoAlphabet::global());
    AlignmentSubstMatrix<double> M(matrixFile,alphabet);
    Sequence *seq1=Sequence::load(file1,alphabet);
    Sequence *seq2=Sequence::load(file2,alphabet);

    int n=cmd.option('n') ? cmd.optParm('n').asInt() : 1;
    for(int i=0 ; i<n ; ++i) {
      Needleman<double> aligner(alphabet,*seq1,*seq2,M,gapOpen,gapExtend);
      AlignmentPath *alignment=local ? aligner.partialAlignment() : aligner.fullAlignment(wantUngapped);
      if(cmd.option('c')) cigarFile<<alignment->getCigarString()<<endl;
      int mismatches, insertions;
      alignment->countMismatches(mismatches,insertions);
      int alignmentLength=alignment->getAlignmentLength();
      int matches=alignmentLength-insertions-mismatches;
      float percentIdentity=int(1000*matches/float(alignmentLength))/10;
      double score=alignment->getScore();
      if(quiet) {
	int firstBegin, firstEnd, secondBegin, secondEnd;
	alignment->getMatchExtent(firstBegin,firstEnd,secondBegin,secondEnd);
	cout<<firstBegin<<"\t"<<firstEnd<<"\t"<<secondBegin<<"\t"<<secondEnd
            <<"\t"<<matches<<"\t"<<score<<endl;
      }
      else if(terseUngapped) {
	int firstBegin, firstEnd, secondBegin, secondEnd;
	alignment->getMatchExtent(firstBegin,firstEnd,secondBegin,secondEnd);
	cout<<firstBegin<<"\t"<<firstEnd<<"\t"<<secondBegin<<"\t"<<secondEnd
	    <<"\t"<<matches<<endl;
      }
      else {
	cout << "\n-----------------------------------------------------------------------\n";
	cout << "Sequence #1: " << file1 << endl;
	cout << "Sequence #2: " << file2 << endl;
	cout << "Percent identity: " << percentIdentity
	     << "%, matches=" << matches << ", mismatches=" << mismatches
	     << ",\ninsertions=" << insertions
	     << ", alignment length=" << alignmentLength 
	     << ", score=" << score
	     << endl;
	cout << "\n<alignment>\n\n" << *alignment << "\n</alignment>\n";
	cout << "-----------------------------------------------------------------------\n";
      }
    }
    return 0;
  }

