/****************************************************************
 ProjectionChecker.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/GffReader.H"
#include "BOOM/ProteinTrans.H"
#include "BOOM/CodonIterator.H"
#include "Labeling.H"
using namespace std;
using namespace BOOM;



class ProjectionChecker {
public:
  ProjectionChecker(const GffTranscript &refTrans,
		    const GffTranscript &altTrans,
		    const String &refSubstrate,
		    const String &altSubstrate,
		    const Labeling &labeling,
		    const Set<String> &nonCanonicalGTs, 
		    const Set<String> &nonCanonicalAGs);
  bool checkSpliceSites(bool quiet);
  bool checkFrameshifts(const Labeling &,const GffTranscript &,
			const String &substrate,bool quiet);
  bool detectNMD(GffTranscript &altTrans,const String &altSubstrate,
		 bool quiet);
  void translate(GffTranscript &refTrans,GffTranscript &altTrans,
		 String &refProtein,String &altProtein);
  bool hasStartCodon(const String &protein);
  bool hasStopCodon(const String &protein);
  bool hasPTC(const String &protein,int &PTCpos);
  static bool geneIsWellFormed(const GffTranscript &,
			       const String &substrate,
			       const Set<String> &noncanonicalGTs,
			       const Set<String> &noncanonicalAGs,
			       bool &noStart,bool &noStop,                                                    
			       bool &PTC,bool &badSpliceSite,
			       String &msg);
protected:
  const Set<String> &nonCanonicalGTs, &nonCanonicalAGs;
  const GffTranscript &refTrans, &altTrans;
  const String &refSubstrate, &altSubstrate;
  const Labeling &labeling;
  bool checkDonor(GffExon &refExon,const String &refSubstrate,
		  GffExon &altExon,const String &altSubstrate,bool quiet);
  bool checkAcceptor(GffExon &refExon,const String &refSubstrate,
		     GffExon &altExon,const String &altSubstrate,bool quiet);
  static bool checkDonor(GffExon &exon,const String &substrate,
		  const Set<String> &noncanoncialGTs,String &msg);
  static bool checkAcceptor(GffExon &exon,const String &substrate,
			    const Set<String> &noncanonicalAGs,String &msg);
  static String getDonor(GffExon &,const String &substrate,int &pos);
  static String getAcceptor(GffExon &,const String &substrate,int &pos);
};



