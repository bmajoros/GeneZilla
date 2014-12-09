/**************************************************************
 Strand.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#ifndef INCL_strand_H
#define INCL_strand_H

typedef char Strand;
inline Strand complement(Strand s) {return s=='-' ? '+' : '-';}

#define PLUS_STRAND '+'
#define MINUS_STRAND '-'
#define FORWARD_STRAND PLUS_STRAND
#define REVERSE_STRAND MINUS_STRAND
#define EITHER_STRAND '.'
#define NO_STRAND '.'

#endif

