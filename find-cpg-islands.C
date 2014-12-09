/****************************************************************
 gilligan.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "gilligan.H"
#include "BOOM/FastaReader.H"
#include "BOOM/IndexedFasta.H"


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
    BOOM::CommandLine cmd(argc,argv,"s:w:i:g:r:c");
    if(cmd.numArgs()!=1)
      throw BOOM::String(
"\ngilligan <*.fasta> [options]\n\
          where:\n\
                -s # = max separation for merging islands (200)\n\
                -w # = sliding window size (100)\n\
                -i # = minimum island size (200)\n\
                -g # = min GC% for an island (0.5)\n\
                -r # = min observed/expected CpG dinucleotide ratio (0.6)\n\
                -c : fasta file is 'compiled' (no defline or whitespace)\n\
");
    maxSeparation = cmd.option('s') ? cmd.optParm('s').asInt() : 200;
    windowSize    = cmd.option('w') ? cmd.optParm('w').asInt() : 100;
    minIslandSize = cmd.option('i') ? cmd.optParm('i').asInt() : 200;
    minGC         = cmd.option('g') ? cmd.optParm('g').asFloat() : 0.5;
    minRatio      = cmd.option('r') ? cmd.optParm('r').asFloat() : 0.6;
    bool isCompiled=cmd.option('c');
    BOOM::String fastaFile=cmd.arg(0);

    if(isCompiled)
      {
      }
    else
      {
	// Iterate through all the sequences in the multi-FASTA file
	BOOM::FastaReader reader(fastaFile);
	BOOM::String dummy;
	while(reader.nextSequence(defline,sequence))
	  {
	    seqLen=sequence.length();
	    if(seqLen<windowSize || seqLen<minIslandSize) continue;
	    seqStr=sequence.c_str();
	    BOOM::FastaReader::parseDefline(defline,substrateID,dummy);
	    
	    // Make a first pass over the sequence using the sliding window
	    // to produce an initial estimate of where the CpG islands are 
	    // located
	    fillArray();
	    
	    // Erase any putative island which is too short
	    eraseSmallIslands();
	    
	    // Merge together any islands which are closer than some 
	    // threshold
	    mergeCloseIslands();
	    
	    // Output the remaining islands in GFF format
	    generateOutput();
	  }
      }

    return 0;
  }



void Application::fillArray()
{
  // Initialize the array
  array.resize(seqLen);
  array.setAllTo(ZERO);
  
  // Initialize variables & counts for first window
  int end=windowSize-1;
  int C=0, G=0, CpG=0;
  char prevBase='?';
  const char *p=seqStr, *pEnd=seqStr+windowSize;
  for(; p!=pEnd ; ++p)
    {
      char x=*p;
      switch(x)
	{
	case 'C': ++C; break;
	case 'G': ++G; if(prevBase=='C') ++CpG; break;
	}
      prevBase=x;
    }

  // Now slide the window along the rest of the sequence
  const float N=windowSize;
  bool prevWindowWasMarked=false;
  const char *pBegin=seqStr;
  --pEnd;
  while(true)
    {
      // First, handle the marking of this window
      float GCpercent=(C+G)/N;
      float ratio=CpG*N/(C*G);
      //cout<<GCpercent<<"\t"<<ratio<<endl;//###
      if(GCpercent>=minGC && ratio>=minRatio)
	if(prevWindowWasMarked)
	  array[end]=ONE;
	else
	  {
	    for(int i=end-windowSize+1 ; i<=end ; ++i)
	      array[i]=ONE;
	    prevWindowWasMarked=true;
	  }
      
      // Now update all variables for the next window
      ++end;
      if(end>=seqLen) break;
      switch(*pBegin)
	{
	case 'C': --C; if(*(pBegin+1)=='G') --CpG; break;
	case 'G': --G; break;
	}
      ++pBegin;
      ++pEnd;
      char x=*pEnd;
      switch(x)
	{
	case 'C': ++C; break;
	case 'G': ++G; if(prevBase=='C') ++CpG; break;
	}
      prevBase=x;
    }
}



void Application::eraseSmallIslands()
{
  int islandLength=0;
  for(int i=0 ; i<seqLen ; ++i)
    if(array[i]==ONE) ++islandLength;
    else
      {
	if(islandLength>0 && islandLength<minIslandSize)
	  for(int j=i-islandLength ; j<i ; ++j)
	    array[j]=ZERO;
	islandLength=0;
      }
  if(islandLength>0 && islandLength<minIslandSize)
    for(int j=seqLen-islandLength ; j<seqLen ; ++j)
      array[j]=ZERO;
}



void Application::mergeCloseIslands()
{
  int seaLength=0;
  for(int i=0 ; i<seqLen ; ++i)
    if(array[i]==ZERO) ++seaLength;
    else
      {
	if(seaLength>0 && seaLength<maxSeparation)
	  for(int j=i-seaLength ; j<i ; ++j)
	    array[j]=ONE;
	seaLength=0;
      }
}



void Application::generateOutput()
{
  int islandLength=0;
  for(int i=0 ; i<seqLen ; ++i)
    if(array[i]==ONE) ++islandLength;
    else
      {
	if(islandLength>0)
	  {
	    int begin=i-islandLength+1;
	    int end=i;
	    cout<<substrateID<<"\t"<<"CpG\tisland\t"<<begin<<"\t"<<end
		<<"\t"<<".\t+\t."<<endl;
	  }
	islandLength=0;
      }
}



