/****************************************************************
 ContentRegions.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "ContentRegions.H"
using namespace std;
using namespace BOOM;

ContentRegions::ContentRegions(const GffTranscript &transcript,int seqLength)
{
  init(transcript,seqLength);
}



const Vector<ContentRegion> &ContentRegions::asVector()
{
  return regions;
}



void ContentRegions::init(const GffTranscript &transcript,int seqLength)
{
  const int numExons=transcript.numExons();
  if(numExons<1) throw "Empty transcript in ContentRegions::init()";
  regions.push_back(ContentRegion(INTERGENIC,0,transcript.getBegin()));
  for(int i=0 ; i<numExons ; ++i) {
    const GffExon &exon=transcript.getIthExon(i);
    ContentType type;
    if(numExons==1) type=SINGLE_EXON;
    else if(i==0) type=INITIAL_EXON;
    else if(i==numExons-1) type=FINAL_EXON;
    else type=INTERNAL_EXON;
    regions.push_back(ContentRegion(type,exon.getBegin(),exon.getEnd()));
    if(i+1<numExons) {
      const GffExon &nextExon=transcript.getIthExon(i+1);
      regions.push_back(ContentRegion(INTRON,exon.getEnd(),nextExon.getBegin()));
    }
  }
  regions.push_back(ContentRegion(INTERGENIC,transcript.getEnd(),seqLength));
}



bool ContentRegions::findJunction(int pos,const ContentRegion *&preceding,
				  const ContentRegion *&following) const
{
  for(Vector<ContentRegion>::const_iterator cur=regions.begin(), 
	end=regions.end() ; cur!=end ; ++cur) {
    const ContentRegion *region=&*cur;
    if(region->getInterval().getEnd()==pos) {
      preceding=region;
      ++cur;
      following=&*cur;
    }
  }
  return false;
}



ContentRegion *ContentRegions::regionOverlapping(int pos) const
{
  for(Vector<ContentRegion>::const_iterator cur=regions.begin(), 
	end=regions.end() ; cur!=end ; ++cur)
    if((*cur).getInterval().contains(pos)) return &*cur;
  return NULL;
}





