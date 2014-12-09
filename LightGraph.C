/****************************************************************
 LightGraph.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <fstream>
#include <iostream>
#include "LightGraph.H"
#include "BOOM/Time.H"
using namespace std;
using namespace BOOM;

LightGraph::LightGraph(const String &filename)
{
  File f(filename);
  load(f);
}



LightGraph::LightGraph(File &f)
{
  load(f);
}



LightGraph::~LightGraph()
{
  Vector<LightVertex*>::iterator vCur=vertices.begin(), vEnd=vertices.end();
  for(; vCur!=vEnd ; ++vCur) delete *vCur;
  Vector<LightEdge*>::iterator eCur=edges.begin(), eEnd=edges.end();
  for(; eCur!=eEnd ; ++eCur) delete *eCur;
}



bool LightGraph::save(const String &filename)
{
  ofstream os(filename.c_str());
  save(os);
}



bool LightGraph::save(ostream &os)
{
  os.precision(16);
  const String timestamp=getDateAndTime();
  const int V=vertices.size(), E=edges.size();
  os<<"##gff-version 2\n\
##source-version reweight-graph ver1.x\n\
##date "<<timestamp<<"\n\
##Type ORFgraph\n\
# stats: "<<V<<" vertices, "<<E<<" edges, "<<substrateLength<<" residues\n\
";
  for(int i=0 ; i<V ; ++i) {
    LightVertex *v=vertices[i];
    if(!v) continue;
    String vertexType;
    if(v->getBegin()<0) vertexType="left_terminus";
    else if(v->getBegin()>=substrateLength) vertexType="right_terminus";
    else vertexType="vertex";
    v->printOn(os,vertexType);
    os<<endl;
  }
  for(int i=0 ; i<E ; ++i) {
    LightEdge *edge=edges[i];
    if(!edge) continue;
    os<<*edge;
  }	
}



int LightGraph::getNumVertices() const
{
  return vertices.size();
}



int LightGraph::getNumEdges() const
{
  return edges.size();
}



LightVertex *LightGraph::getVertex(int i)
{
  return vertices[i];
}



LightEdge *LightGraph::getEdge(int i)
{
  return edges[i];
}



const String &LightGraph::getSubstrate() const
{
  return substrate;
}



void LightGraph::printOn(ostream &os)
{
  save(os);
}



bool LightGraph::load(File &f)
{
  bool seenLeft=false, seenRight=false;
  int lastEdgeID=-1;
  while(!f.eof()) {
    String line=f.getline();
    line.trimWhitespace();
    BOOM::Vector<BOOM::String> &fields=*line.getFields();
    int numFields=fields.size();
    bool keep=true;
    if(numFields>=8) {
      const String recType=fields[1];
      if(recType=="stats:") substrateLength=fields[6].asInt();
      else {
	substrate=fields[0];
	if(recType=="vertex" || recType=="left_terminus" ||
	   recType=="right_terminus") {

	  if(recType=="left_terminus") {
	    if(seenLeft) keep=false;
	    seenLeft=true;
	  }
	  else if(recType=="right_terminus") {
	    if(seenRight) keep=false;
	    seenRight=true;
	  }

	  SignalType sigType=stringToSignalType(fields[2]);
	  int begin=fields[3].asInt()-1;
	  int end=fields[4].asInt();
	  float score=fields[5].asFloat();
	  char strand=fields[6][0];
	  if(strand=='-') sigType=reverseComplement(sigType);
	  int ID=vertices.size();
	  LightVertex *v=
	    keep?
	    new LightVertex(substrate,sigType,begin,end,score,strand,ID)
	    : NULL;
	  if(v) {
	    bool support=false;
	    if(fields.size()>=12) support=fields[11].substring(4).asInt();
	    v->setSupport(support);
	  }
	  vertices.push_back(v);
	}
	else if(recType=="edge") {
	  ContentType conType=stringToContentType(fields[2]);
	  int begin=fields[3].asInt()-1;
	  int end=fields[4].asInt();
	  float score=fields[5].asFloat();
	  char strand=fields[6][0];
	  if(strand=='-') conType=reverseComplement(conType);
	  int frame=fields[7].asInt();
	  int leftID=fields[8].substring(5).asInt();
	  int rightID=fields[9].substring(6).asInt();
	  int edgeID=fields[10].substring(7).asInt();
	  bool support=false;
	  if(fields.size()>=12) support=fields[11].substring(4).asInt();
	  LightVertex *left=vertices[leftID], *right=vertices[rightID];
	  if(!left || !right) { delete &fields; continue; }
	  if(edgeID!=lastEdgeID) {
	    LightEdge *edge=
	      left && right ? 
	      new LightEdge(substrate,conType,left,right,begin,end,
			    strand,edgeID)
	      : NULL;
	    if(edge) edge->setSupport(support);
	    edges.push_back(edge);
	    if(left) left->addEdgeOut(edge);
	    if(right) right->addEdgeIn(edge);
	    lastEdgeID=edgeID;
	  }
	  edges[edges.size()-1]->setScore(frame,score);
	}	
      }
    }
    delete &fields;
  }
}



ostream &operator<<(ostream &os,const LightGraph &g)
{
  g.printOn(os);
  return os;
}



int LightGraph::getLargestEdgeID() const
{
  int largest=-1;

  for(Vector<LightEdge*>::const_iterator cur=edges.begin(), end=edges.end() ; 
      cur!=end ; ++cur) {
    int id=(*cur)->getID();
    if(id>largest) largest=id;
  }
  return largest;
}


