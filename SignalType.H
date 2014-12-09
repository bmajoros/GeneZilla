/*
 SignalType.*
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#ifndef INCL_SignalType_H
#define INCL_SignalType_H

#include <iostream>
#include "BOOM/String.H"
#include "Strand.H"

enum SignalType
  {
    ATG,
    TAG,
    GT,
    AG,
    GT_U12,
    AG_U12,
    PROM,
    POLYA,
    NEG_ATG,
    NEG_TAG,
    NEG_GT,
    NEG_AG,
    NEG_GT_U12,
    NEG_AG_U12,
    NEG_PROM,
    NEG_POLYA,

    // Extended types (used only for transition probabilities):
    GT0,
    GT1,
    GT2,
    AG0,
    AG1,
    AG2,
    NEG_GT0,
    NEG_GT1,
    NEG_GT2,
    NEG_AG0,
    NEG_AG1,
    NEG_AG2,

    NO_SIGNAL_TYPE
  };

inline int numSignalTypes() {return 12;} // doesn't include extended types!
inline SignalType &operator++(SignalType &s) {s=SignalType(int(s)+1);return s;}

SignalType dropStrand(SignalType);
ostream &operator<<(ostream &,SignalType);
istream &operator>>(istream &,SignalType &);
BOOM::String signalTypeToString(SignalType);
BOOM::String signalTypeToName(SignalType);
SignalType reverseComplement(SignalType);
SignalType stringToSignalType(const BOOM::String &);
Strand getStrand(SignalType);
bool endsCoding(SignalType);
bool beginsCoding(SignalType);
bool beginsIntron(SignalType);
bool endsIntron(SignalType);

class TypeNamesInitializer
{
  static TypeNamesInitializer initializer;
  TypeNamesInitializer();
};

#endif
