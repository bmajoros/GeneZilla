/*
 ContentType.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 */
#include "ContentType.H"
#include "BOOM/Map.H"


static BOOM::Map<BOOM::String,ContentType> contentTypeMap;
static BOOM::Map<BOOM::String,ContentType> longContentTypeMap;
ContentTypeInitializer ContentTypeInitializer::initializer;


ContentTypeInitializer::ContentTypeInitializer()
{
  contentTypeMap["INITIAL-EXON"]=INITIAL_EXON;
  contentTypeMap["INTERNAL-EXON"]=INTERNAL_EXON;
  contentTypeMap["FINAL-EXON"]=FINAL_EXON;
  contentTypeMap["SINGLE-EXON"]=SINGLE_EXON;
  contentTypeMap["INTRON"]=INTRON;
  contentTypeMap["INTERGENIC"]=INTERGENIC;
  contentTypeMap["FIVE-PRIME-UTR"]=FIVE_PRIME_UTR;
  contentTypeMap["THREE-PRIME-UTR"]=THREE_PRIME_UTR;
  contentTypeMap["NEG-INITIAL-EXON"]=NEG_INITIAL_EXON;
  contentTypeMap["NEG-INTERNAL-EXON"]=NEG_INTERNAL_EXON;
  contentTypeMap["NEG-FINAL-EXON"]=NEG_FINAL_EXON;
  contentTypeMap["NEG-SINGLE-EXON"]=NEG_SINGLE_EXON;
  contentTypeMap["NEG-INTRON"]=NEG_INTRON;
  contentTypeMap["NEG-FIVE-PRIME-UTR"]=NEG_FIVE_PRIME_UTR;
  contentTypeMap["NEG-THREE-PRIME-UTR"]=NEG_THREE_PRIME_UTR;
  contentTypeMap["UNKNOWN-CONTENT-FORWARD"]=UNKNOWN_CONTENT_FORWARD;
  contentTypeMap["UNKNOWN-CONTENT-REVERSE"]=UNKNOWN_CONTENT_REVERSE;

  longContentTypeMap["initial-exon"]=INITIAL_EXON;
  longContentTypeMap["internal-exon"]=INTERNAL_EXON;
  longContentTypeMap["final-exon"]=FINAL_EXON;
  longContentTypeMap["single-exon"]=SINGLE_EXON;
  longContentTypeMap["intron"]=INTRON;
  longContentTypeMap["intergenic"]=INTERGENIC;
  longContentTypeMap["5'-UTR"]=FIVE_PRIME_UTR;
  longContentTypeMap["3'-UTR"]=THREE_PRIME_UTR;
}



Strand getStrand(ContentType t)
{
  switch(t)
    {
    case UNKNOWN_CONTENT_FORWARD:
    case INITIAL_EXON:
    case INTERNAL_EXON:
    case FINAL_EXON:
    case SINGLE_EXON:
    case INTRON:
    case INTERGENIC://###
    case FIVE_PRIME_UTR:
    case THREE_PRIME_UTR:
      return FORWARD_STRAND;

    case UNKNOWN_CONTENT_REVERSE:
    case NEG_INITIAL_EXON:
    case NEG_INTERNAL_EXON:
    case NEG_FINAL_EXON:
    case NEG_SINGLE_EXON:
    case NEG_INTRON:
    case NEG_FIVE_PRIME_UTR:
    case NEG_THREE_PRIME_UTR:
      return REVERSE_STRAND;

    default: throw BOOM::String(__FILE__)+__LINE__;
    }
}



BOOM::String contentTypeToString(ContentType t)
{
  switch(t)
    {
    case UNKNOWN_CONTENT_FORWARD: return "UNKNOWN-CONTENT-FORWARD";
    case UNKNOWN_CONTENT_REVERSE: return "UNKNOWN-CONTENT-REVERSE";
    case INITIAL_EXON:         return "INITIAL-EXON";
    case INTERNAL_EXON:        return "INTERNAL-EXON";
    case FINAL_EXON:           return "FINAL-EXON";
    case SINGLE_EXON:          return "SINGLE-EXON";
    case INTRON:               return "INTRON";
    case INTERGENIC:           return "INTERGENIC";
    case FIVE_PRIME_UTR:       return "FIVE-PRIME-UTR";
    case THREE_PRIME_UTR:      return "THREE-PRIME-UTR";
    case NEG_INITIAL_EXON:     return "NEG-INITIAL-EXON";
    case NEG_INTERNAL_EXON:    return "NEG-INTERNAL-EXON";
    case NEG_FINAL_EXON:       return "NEG-FINAL-EXON";
    case NEG_SINGLE_EXON:      return "NEG-SINGLE-EXON";
    case NEG_INTRON:           return "NEG-INTRON";
    case NEG_FIVE_PRIME_UTR:   return "NEG-FIVE-PRIME-UTR";
    case NEG_THREE_PRIME_UTR:  return "NEG-THREE-PRIME-UTR";

    default: throw BOOM::String("INVALID CONTENT TYPE: ")+int(t);
    }
}



BOOM::String contentTypeNiceString(ContentType t)
{
  switch(t)
    {
    case UNKNOWN_CONTENT_FORWARD: return "UNKNOWN-CONTENT-FORWARD";
    case UNKNOWN_CONTENT_REVERSE: return "UNKNOWN-CONTENT-REVERSE";
    case INITIAL_EXON:         return "initial-exon";
    case INTERNAL_EXON:        return "internal-exon";
    case FINAL_EXON:           return "final-exon";
    case SINGLE_EXON:          return "single-exon";
    case INTRON:               return "intron";
    case INTERGENIC:           return "intergenic";
    case FIVE_PRIME_UTR:       return "5'-UTR";
    case THREE_PRIME_UTR:      return "3'-UTR";
    case NEG_INITIAL_EXON:     return "initial-exon";
    case NEG_INTERNAL_EXON:    return "internal-exon";
    case NEG_FINAL_EXON:       return "final-exon";
    case NEG_SINGLE_EXON:      return "single-exon";
    case NEG_INTRON:           return "intron";
    case NEG_FIVE_PRIME_UTR:   return "5'-UTR";
    case NEG_THREE_PRIME_UTR:  return "3'-UTR";

    default: throw BOOM::String("INVALID CONTENT TYPE: ")+int(t);
    }
}



ContentType stringToContentType(const BOOM::String &str)
{
  if(contentTypeMap.isDefined(str)) return contentTypeMap[str];
  if(longContentTypeMap.isDefined(str)) return longContentTypeMap[str];
  throw BOOM::String("undefined ContentType: ")+str;
}



istream &operator>>(istream &is,ContentType &t)
{
  BOOM::String buf;
  is >> buf;
  t=stringToContentType(buf);
  return is;
}



ostream &operator<<(ostream &os,const ContentType &t)
{
  os << contentTypeToString(t);
  return os;
}



ContentType reverseComplement(ContentType t)
{
  switch(t)
    {
    case INITIAL_EXON:         return NEG_INITIAL_EXON;
    case INTERNAL_EXON:        return NEG_INTERNAL_EXON;
    case FINAL_EXON:           return NEG_FINAL_EXON;
    case SINGLE_EXON:          return NEG_SINGLE_EXON;
    case INTRON:               return NEG_INTRON;
    case FIVE_PRIME_UTR:       return NEG_FIVE_PRIME_UTR;
    case THREE_PRIME_UTR:      return NEG_THREE_PRIME_UTR;
    case UNKNOWN_CONTENT_FORWARD: return UNKNOWN_CONTENT_REVERSE;

    case NEG_INITIAL_EXON:     return INITIAL_EXON;
    case NEG_INTERNAL_EXON:    return INTERNAL_EXON;
    case NEG_FINAL_EXON:       return FINAL_EXON;
    case NEG_SINGLE_EXON:      return SINGLE_EXON;
    case NEG_INTRON:           return INTRON;
    case NEG_FIVE_PRIME_UTR:   return FIVE_PRIME_UTR;
    case NEG_THREE_PRIME_UTR:  return THREE_PRIME_UTR;
    case UNKNOWN_CONTENT_REVERSE: return UNKNOWN_CONTENT_FORWARD;

    case INTERGENIC:           return INTERGENIC;

    default: throw BOOM::String(__FILE__)+__LINE__;
    }
}



bool isCoding(ContentType t)
{
  switch(t)
    {
    case INITIAL_EXON:         
    case INTERNAL_EXON:        
    case FINAL_EXON:           
    case SINGLE_EXON:          
    case NEG_INITIAL_EXON:     
    case NEG_INTERNAL_EXON:    
    case NEG_FINAL_EXON:       
    case NEG_SINGLE_EXON:      
      return true;

    case INTRON:               
    case INTERGENIC:           
    case FIVE_PRIME_UTR:       
    case THREE_PRIME_UTR:      
    case NEG_INTRON:           
    case NEG_FIVE_PRIME_UTR:   
    case NEG_THREE_PRIME_UTR:  
      return false;

    case UNKNOWN_CONTENT_FORWARD:
    case UNKNOWN_CONTENT_REVERSE:
      throw "isCoding(UNKNOWN_CONTENT_FORWARD/REVERSE)";

    default: throw BOOM::String(int(t))+" "+__FILE__+" "+__LINE__;
    }
}



bool isIntron(ContentType t)
{
  switch(t)
    {
    case INTRON:               
    case NEG_INTRON:           
      return true;
    default:
      return false;
    }
}



bool isIntergenic(ContentType t)
{
  return t==INTERGENIC;
}



bool isUTR(ContentType t)
{
  switch(t)
    {
    case FIVE_PRIME_UTR:       
    case THREE_PRIME_UTR:      
    case NEG_FIVE_PRIME_UTR:   
    case NEG_THREE_PRIME_UTR:  
      return true;

    default: 
      return false;
    }
}



SignalType leftSignal(ContentType contentType)
{
  switch(contentType)
    {
    case INITIAL_EXON:             return ATG;
    case INTERNAL_EXON:            return AG;
    case FINAL_EXON:               return AG;
    case SINGLE_EXON:              return ATG;
    case INTRON:                   return GT;
    case FIVE_PRIME_UTR:           return PROM;
    case THREE_PRIME_UTR:          return TAG;
    case NEG_INITIAL_EXON:         return NEG_GT;
    case NEG_INTERNAL_EXON:        return NEG_GT;
    case NEG_FINAL_EXON:           return NEG_TAG;
    case NEG_SINGLE_EXON:          return NEG_TAG;
    case NEG_INTRON:               return NEG_AG;
    case NEG_FIVE_PRIME_UTR:       return NEG_ATG;
    case NEG_THREE_PRIME_UTR:      return NEG_POLYA;
    default: throw BOOM::String(__FILE__)+__LINE__;
    }
}



SignalType rightSignal(ContentType contentType)
{
  switch(contentType)
    {
    case INITIAL_EXON:             return GT;
    case INTERNAL_EXON:            return GT;
    case FINAL_EXON:               return TAG;
    case SINGLE_EXON:              return TAG;
    case INTRON:                   return AG;
    case FIVE_PRIME_UTR:           return ATG;
    case THREE_PRIME_UTR:          return POLYA;
    case NEG_INITIAL_EXON:         return NEG_ATG;
    case NEG_INTERNAL_EXON:        return NEG_AG;
    case NEG_FINAL_EXON:           return NEG_AG;
    case NEG_SINGLE_EXON:          return NEG_ATG;
    case NEG_INTRON:               return NEG_GT;
    case NEG_FIVE_PRIME_UTR:       return NEG_PROM;
    case NEG_THREE_PRIME_UTR:      return NEG_TAG;
    default: throw BOOM::String(__FILE__)+__LINE__;
    }
}



ContentType exonTypeToContentType(ExonType exonType,
				  Strand strand)
{
  switch(exonType)
    {
    case ET_INITIAL_EXON:
      return strand==FORWARD_STRAND ? INITIAL_EXON : NEG_INITIAL_EXON;
    case ET_INTERNAL_EXON:
      return strand==FORWARD_STRAND ? INTERNAL_EXON : NEG_INTERNAL_EXON;
    case ET_FINAL_EXON:
      return strand==FORWARD_STRAND ? FINAL_EXON : NEG_FINAL_EXON;
    case ET_SINGLE_EXON:
      return strand==FORWARD_STRAND ? SINGLE_EXON : NEG_SINGLE_EXON;
    case ET_EXON:       
    default:
      throw "::exonTypeToContentType()";
    }
}



ContentType dropStrand(ContentType t)
{
  switch(t)
    {
    case INITIAL_EXON:
    case INTERNAL_EXON:
    case FINAL_EXON:
    case SINGLE_EXON:
    case INTRON:
    case FIVE_PRIME_UTR:
    case THREE_PRIME_UTR:
    case UNKNOWN_CONTENT_FORWARD:
      return t;
    case NEG_INITIAL_EXON:
    case NEG_INTERNAL_EXON:
    case NEG_FINAL_EXON:
    case NEG_SINGLE_EXON:
    case NEG_INTRON:
    case NEG_FIVE_PRIME_UTR:
    case NEG_THREE_PRIME_UTR:
    case UNKNOWN_CONTENT_REVERSE:
      return reverseComplement(t);
    case INTERGENIC:           
      return INTERGENIC;
    }
  INTERNAL_ERROR;
}


