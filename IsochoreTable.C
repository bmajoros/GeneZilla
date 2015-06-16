/****************************************************************
 IsochoreTable.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "IsochoreTable.H"
#include "BOOM/File.H"
#include "SignalTypeProperties.H"
#include "TopologyLoader.H"
#include "EmpiricalDistribution.H"
#include "GeometricDistribution.H"
#include <iostream>
#include <fstream>
using namespace std;



BOOM::String IsochoreTable::romanNumerals[]=
  {"0","I","II","III","IV","V","VI","VII","VIII","IX","X",
   "XI","XII","XIII","XIV","XV","XVI","XVII","XVIII","XIX",
   "XX"};


// G+C <= 100% : genezilla.cfg
BOOM::Regex IsochoreTable::lineRegex(
  "G[+]C\\s*<=\\s*(\\S+)%\\s*:\\s*(\\S+)");

// GLOBAL three-prime-optimism = 6
BOOM::Regex IsochoreTable::globalRegex(
  "GLOBAL\\s+(\\S+)\\s*=\\s*(\\S+)");

// # this is a comment
BOOM::Regex IsochoreTable::commentRegex("^\\s*#");


IsochoreTable::IsochoreTable(GarbageCollector &garbageCollector)
  : garbageCollector(garbageCollector)
{
  // ctor
}



void IsochoreTable::load(const BOOM::String &filename)
{
  // This method loads the *.iso file and all the config files (*.cfg)
  // that are named in that file, as well as any GLOBAL attributes
  // set in the *.iso file (which supercede those set in the individual
  // *.cfg files).

  // First, parse the *.iso file
  int isochoreId=1;
  BOOM::File infile(filename);
  BOOM::Vector< pair<BOOM::String,BOOM::String> > pairs;
  BOOM::Vector<BOOM::ConfigFile*> configFiles;
  while(!infile.eof())
    {
      BOOM::String line=infile.readLine();
      line.trimWhitespace();
      if(lineRegex.search(line))
	{
	  double maximumGC=lineRegex[1].asDouble()/100;
	  BOOM::String configFilename=lineRegex[2];
	  BOOM::ConfigFile *configFile=new ConfigFile(configFilename);
	  configFiles.push_back(configFile);
	  Isochore *isochore=new Isochore;
	  isochore->configFile=*configFile;
	  isochore->name=toRomanNumerals(isochoreId);
	  isochore->maximumGC=maximumGC;
	  isochores.push_back(isochore);
	  boundaries.push_back(maximumGC);
	  //loadIsochore(configFile,*isochore);
	  ++isochoreId;
	}
      else if(globalRegex.search(line))
	{
	  BOOM::String key=globalRegex[1];
	  BOOM::String value=globalRegex[2];
	  pairs.push_back(pair<BOOM::String,BOOM::String>(key,value));
	}
      else if(commentRegex.search(line)) continue;
      else if(!line.isEmpty())
	throw BOOM::String("Syntax error in ISO file:\n")+line;
    }
  infile.close();

  // Add the GLOBAL attributes to all the config file representations
  int n=pairs.size();
  int m=isochores.size();
  for(int i=0 ; i<n ; ++i)
    {
      pair<BOOM::String,BOOM::String> &p=pairs[i];
      for(int j=0 ; j<m ; ++j)
	configFiles[j]->enter(p.first,p.second);
      for(Vector<Isochore*>::iterator cur=isochores.begin(), end=isochores.end() ;
	  cur!=end ; ++cur)
	(*cur)->configFile.enter(p.first,p.second);
    }

  // Finally, extract the gene-finder parameters from the config files
  for(int j=0 ; j<m ; ++j)
    loadIsochore(*configFiles[j],*isochores[j]);
}



BOOM::String IsochoreTable::toRomanNumerals(int x)
{
  if(x<sizeof(romanNumerals)/sizeof(romanNumerals[0]))
    return romanNumerals[x];
  return BOOM::String("")+x;
}



Isochore *IsochoreTable::getIsochore(float GC)
{
  int n=boundaries.size();
  for(int i=0 ; i<n ; ++i)
    if(GC<boundaries[i])
      return isochores[i];
  for(int i=0 ; i<n ; ++i) cout<<boundaries[i]<<endl;
  throw String("GC% ")+GC+" not found in IsochoreTable::getIsochore()";
}



Isochore *IsochoreTable::getIthIsochore(int i)
{
  return isochores[i];
}



int IsochoreTable::getQueueCapacity()
{
  return queueCapacity;
}



void IsochoreTable::loadIsochore(BOOM::ConfigFile &configFile,Isochore &iso)
{
  iso.threePrimeOptimism=
    configFile.lookup("three-prime-optimism").asDouble();
  iso.fivePrimeOptimism=
    configFile.lookup("five-prime-optimism").asDouble();
  iso.interpolateHistograms=
    configFile.getBoolOrDie("histogram-interpolation");
  queueCapacity=configFile.getIntOrDie("queue-capacity");
  iso.useSignalThresholds=configFile.getBoolOrDie("use-signal-thresholds");
  iso.signalThresholdMultiplier=
    configFile.getDoubleOrDie("signal-threshold-multiplier");
  BOOM::String topologyFile=configFile.lookupOrDie("metamodel-topology");
  if(!SignalTypeProperties::global.hasBeenLoaded())
    TopologyLoader::load(topologyFile);
  iso.pI=configFile.getDoubleOrDie("probability-start-in-intron");
  iso.pN=configFile.getDoubleOrDie("probability-start-in-intergenic");
  iso.pF=configFile.getDoubleOrDie("probability-start-in-5pUTR");
  iso.pT=configFile.getDoubleOrDie("probability-start-in-3pUTR");
  double pSum=iso.pI+iso.pN+iso.pF+iso.pT;
  iso.pI/=pSum; iso.pN/=pSum; iso.pF/=pSum; iso.pT/=pSum;
  
  loadSubmodels(configFile,iso);

  BOOM::String transFile=configFile.lookupOrDie("transition-probabilities");
  float optimism=configFile.getFloatOrDie("exon-optimism");
  float intronOptimism=configFile.getFloatOrDie("intron-optimism");
  ifstream is(transFile.c_str());
  if(!is.good()) throw BOOM::String("Can't open file ")+transFile;
  int numSignalTypes;
  is >> numSignalTypes;
  iso.transitionProbs=
    new Transitions(numSignalTypes,is,optimism,intronOptimism);
}



void IsochoreTable::loadSubmodels(BOOM::ConfigFile &configFile,
				  Isochore &iso)
{
  // Load signal sensors
  loadSignalSensor(configFile,"donor-model","donor-consensus",iso);
  if(configFile.isDefined("donor-model-U12"))
     loadSignalSensor(configFile,"donor-model-U12","donor-consensus-U12",iso);
  loadSignalSensor(configFile,"acceptor-model","acceptor-consensus",iso);
  if(configFile.isDefined("acceptor-model-U12"))
     loadSignalSensor(configFile,"acceptor-model-U12","acceptor-consensus-U12",
		      iso);
  loadSignalSensor(configFile,"start-codon-model","start-codon-consensus",
		   iso);
  loadSignalSensor(configFile,"stop-codon-model","stop-codon-consensus",
		   iso);
  loadSignalSensor(configFile,"polya-model","polya-consensus",iso);
  loadSignalSensor(configFile,"promoter-model","promoter-consensus",iso);
  BOOM::String shortSignalPeptideFile=
    configFile.lookup("short-signal-peptide-model");
  if(shortSignalPeptideFile!="") 
    loadSignalPeptideModel(shortSignalPeptideFile,configFile,iso);
  BOOM::String longSignalPeptideFile=
    configFile.lookup("long-signal-peptide-model");
  if(longSignalPeptideFile!="") 
    loadSignalPeptideModel(longSignalPeptideFile,configFile,iso);

  // Load content sensors
  loadContentSensor(configFile,"initial-exons",EMPIRICAL_DISTRIBUTION,
		    "initial-exon-lengths",INITIAL_EXON,iso);
  loadContentSensor(configFile,"internal-exons",EMPIRICAL_DISTRIBUTION,
		    "internal-exon-lengths",INTERNAL_EXON,iso);
  loadContentSensor(configFile,"final-exons",EMPIRICAL_DISTRIBUTION,
		    "final-exon-lengths",FINAL_EXON,iso);
  loadContentSensor(configFile,"single-exons",EMPIRICAL_DISTRIBUTION,
		    "single-exon-lengths",SINGLE_EXON,iso);
  loadContentSensor(configFile,"introns",GEOMETRIC_DISTRIBUTION,
		    "mean-intron-length",INTRON,iso);
  loadContentSensor(configFile,"intergenic",GEOMETRIC_DISTRIBUTION,
		    "mean-intergenic-length",INTERGENIC,iso);
  loadContentSensor(configFile,"3'-UTR",GEOMETRIC_DISTRIBUTION,
		    "mean-3'-UTR-length",THREE_PRIME_UTR,iso);
  loadContentSensor(configFile,"5'-UTR",GEOMETRIC_DISTRIBUTION,
		    "mean-5'-UTR-length",FIVE_PRIME_UTR,iso);
}



void IsochoreTable::loadContentSensor(BOOM::ConfigFile &configFile,
				      const BOOM::String &modelLabel,
				      DistributionType distrType,
				      const BOOM::String &distrLabel,
				      ContentType contentType,
				      Isochore &iso)
{
  // Load the model
  ContentType revContentType=::reverseComplement(contentType);
  BOOM::String filename=configFile.lookupOrDie(modelLabel);
  ContentSensor *sensor;
  if(iso.uniqContentSensors.isDefined(filename)) 
    sensor=iso.uniqContentSensors[filename];
  else 
    {
      sensor=ContentSensor::load(filename);
      iso.uniqContentSensors[filename]=sensor;
      iso.contentSensors.push_back(sensor);

      ContentSensor *rev=sensor->reverseComplement();
      if(rev && rev!=sensor) iso.contentSensors.push_back(rev);
    }
  iso.contentToSensor[contentType]=sensor;
  ContentSensor *rev=sensor->reverseComplement();
  if(rev) iso.contentToSensor[revContentType]=rev;

  // Load or create the distribution
  DiscreteDistribution *distribution;
  BOOM::String distrParm=configFile.lookupOrDie(distrLabel);
  if(distrType==EMPIRICAL_DISTRIBUTION)
    {
      if(!iso.uniqDistributions.isDefined(distrParm))
	iso.uniqDistributions[distrParm]=
	  new EmpiricalDistribution(distrParm,iso.interpolateHistograms);
      distribution=iso.uniqDistributions[distrParm];
    }
  else
    distribution=new GeometricDistribution(distrParm.asInt());
  iso.contentToDistribution[contentType]=distribution;
  iso.contentToDistribution[revContentType]=distribution;

  // Populate the associative arrays of signal comparators for the various
  // queue types:
  if(::isIntron(contentType))
    {
      BOOM::Array1D<SinglePhaseComparator*> &array=*
	new BOOM::Array1D<SinglePhaseComparator*>(3);
      for(int i=0 ; i<3 ; ++i)
	array[i]=new SinglePhaseComparator(i,contentType,*distribution);
      iso.intronComparators[contentType]=&array;

      BOOM::Array1D<SinglePhaseComparator*> &revArray=*
	new BOOM::Array1D<SinglePhaseComparator*>(3);
      for(int i=0 ; i<3 ; ++i)
	revArray[i]=
	  new SinglePhaseComparator(i,revContentType,*distribution);
      iso.intronComparators[revContentType]=&revArray;
    }
  else if(!::isCoding(contentType))
    {
      iso.noncodingComparators[contentType]=
	new NoncodingComparator(contentType,*distribution);
      iso.noncodingComparators[revContentType]=
	new NoncodingComparator(revContentType,*distribution);
    }
}



void IsochoreTable::loadSignalSensor(BOOM::ConfigFile &configFile,
				     const BOOM::String &modelLabel,
				     const BOOM::String &consensusLabel,
				     Isochore &iso)
{
  // Load the model
  BOOM::String filename=configFile.lookupOrDie(modelLabel);
  SignalSensor *sensor=SignalSensor::load(filename,garbageCollector);
  if(!iso.useSignalThresholds) sensor->ignoreCutoff();
  else sensor->multiplyCutoffBy(iso.signalThresholdMultiplier);

  // Load consensus strings
  BOOM::String consensuses=configFile.lookupOrDie(consensusLabel);
  BOOM::Vector<BOOM::String> *fields=consensuses.getFields(",|;:/.&_-+");
  BOOM::Vector<BOOM::String>::iterator cur=fields->begin(), 
    end=fields->end();
  for(; cur!=end ; ++cur)
    sensor->addConsensus(*cur);
  delete fields;

  // Reverse-complement the model
#ifndef FORWARD_STRAND_ONLY
  SignalSensor *negSensor=sensor->reverseComplement();
#endif

  // Add model to global model list
  iso.signalSensors.push_back(sensor);
#ifndef FORWARD_STRAND_ONLY
  iso.signalSensors.push_back(negSensor);
#endif

  // Keep a special pointer to the stop codon sensor (since we
  // will be using it a lot to detect ends of ORFs)
  SignalType signalType=sensor->getSignalType();
  if(signalType==TAG)
    {
      iso.stopCodonSensor=sensor;
#ifndef FORWARD_STRAND_ONLY
      iso.negStopCodonSensor=negSensor;
#endif
    }

  // Set up a mapping from signal types to their sensors
  iso.signalTypeToSensor[signalType]=sensor;
#ifndef FORWARD_STRAND_ONLY
  iso.signalTypeToSensor[negSensor->getSignalType()]=negSensor;
#endif
}



int IsochoreTable::getNumIsochores()
{
  return isochores.size();
}



void IsochoreTable::loadSignalPeptideModel(const BOOM::String &filename,
					   BOOM::ConfigFile &configFile,
					   Isochore &iso)
{
  // Load the model
  SignalSensor *sensor=SignalSensor::load(filename,garbageCollector);
  if(!iso.useSignalThresholds) sensor->ignoreCutoff();
  else sensor->multiplyCutoffBy(iso.signalThresholdMultiplier);

  // Load consensus strings
  BOOM::String consensuses=configFile.lookupOrDie("start-codon-consensus");
  BOOM::Vector<BOOM::String> *fields=consensuses.getFields(",|;:/.&_-+");
  BOOM::Vector<BOOM::String>::iterator cur=fields->begin(), end=fields->end();
  for(; cur!=end ; ++cur)
    sensor->addConsensus(*cur);
  delete fields;

  // Reverse-complement the model
#ifndef FORWARD_STRAND_ONLY
  SignalSensor *negSensor=sensor->reverseComplement();
#endif

  // Add model to global model list
  iso.signalSensors.push_back(sensor);
#ifndef FORWARD_STRAND_ONLY
  iso.signalSensors.push_back(negSensor);
#endif
}



