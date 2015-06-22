/*
 ContentType.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 */
#ifndef INCL_ContentType_H
#define INCL_ContentType_H

#include <iostream>
#include "BOOM/String.H"
#include "BOOM/GffExon.H"
#include "Strand.H"
#include "SignalType.H"
using namespace BOOM;

enum ContentType
{
  UNKNOWN_CONTENT_FORWARD,
  UNKNOWN_CONTENT_REVERSE,
  INITIAL_EXON,
  INTERNAL_EXON,
  FINAL_EXON,
  SINGLE_EXON,
  INTRON,
  INTERGENIC,
  FIVE_PRIME_UTR,
  THREE_PRIME_UTR,
  NEG_INITIAL_EXON,
  NEG_INTERNAL_EXON,
  NEG_FINAL_EXON,
  NEG_SINGLE_EXON,
  NEG_INTRON,
  NEG_FIVE_PRIME_UTR,
  NEG_THREE_PRIME_UTR
};
static const int NumContentTypes=17;

ContentType reverseComplement(ContentType);
ContentType dropStrand(ContentType);
istream &operator>>(istream &,ContentType &);
ostream &operator<<(ostream &,const ContentType &);
Strand getStrand(ContentType);
BOOM::String contentTypeToString(ContentType);
BOOM::String contentTypeNiceString(ContentType);
ContentType stringToContentType(const BOOM::String &);
bool isCoding(ContentType);
bool isIntron(ContentType);
bool isIntergenic(ContentType);
bool isUTR(ContentType);
SignalType leftSignal(ContentType);
SignalType rightSignal(ContentType);
ContentType exonTypeToContentType(ExonType,Strand);
ExonType contentTypeToExonType(ContentType);

class ContentTypeInitializer
{
  static ContentTypeInitializer initializer;
  ContentTypeInitializer();
};

#endif
