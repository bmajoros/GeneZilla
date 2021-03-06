# Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
CC		= g++
DEBUG		= -g
OPTIMIZE	= -O3 -finline-functions -funroll-loops
CFLAGS		= $(OPTIMIZE) -w -fpermissive
LDFLAGS		= $(OPTIMIZE) -w -fpermissive
BOOM		= BOOM
OBJ		= obj
STATIC		=
LIBDIRS		= 
LIBS		= -lpthread  -lgsl -lm -lgslcblas -LBOOM -lBOOM 


all: \
	BOOM \
	obj \
	genezilla \
	crf \
	cia \
	reweight-graph \
	n-best \
	get-bias-model \
	train-signal-sensor \
	train-content-sensor \
	train-branch-acceptor \
	train-tata-cap-model \
	train-signal-peptide-model \
	compile-markov-chain \
	get-transcripts \
	mdd

BOOM:
	ln -s ../BOOM

obj: 
	mkdir obj

.PHONY : clean
clean:
	@rm -f $(OBJ)/*.o

genezilla.H:
	ln -s genezilla-graph.H genezilla.H

#---------------------------------------------------------
$(OBJ)/GffPathFromParseGraph.o:\
		genezilla.H \
		GffPathFromParseGraph.H \
		GffPathFromParseGraph.C
	$(CC) $(CFLAGS) -o $(OBJ)/GffPathFromParseGraph.o -c \
		GffPathFromParseGraph.C
#---------------------------------------------------------
$(OBJ)/RnaJunction.o:\
		RnaJunction.H \
		RnaJunction.C
	$(CC) $(CFLAGS) -o $(OBJ)/RnaJunction.o -c \
		RnaJunction.C
#---------------------------------------------------------
$(OBJ)/RnaJunctions.o:\
		genezilla.H \
		RnaJunctions.H \
		RnaJunctions.C
	$(CC) $(CFLAGS) -o $(OBJ)/RnaJunctions.o -c \
		RnaJunctions.C
#---------------------------------------------------------
$(OBJ)/Propagator.o:\
		genezilla.H \
		Propagator.H \
		Propagator.C
	$(CC) $(CFLAGS) -o $(OBJ)/Propagator.o -c \
		Propagator.C
#---------------------------------------------------------
$(OBJ)/IsochoreTable.o:\
		IsochoreTable.H \
		IsochoreTable.C
	$(CC) $(CFLAGS) -o $(OBJ)/IsochoreTable.o -c \
		IsochoreTable.C
#---------------------------------------------------------
$(OBJ)/CodonTree.o:\
		CodonTree.H \
		CodonTree.C
	$(CC) $(CFLAGS) -o $(OBJ)/CodonTree.o -c \
		CodonTree.C
#---------------------------------------------------------
$(OBJ)/SignalPeptideSensor.o:\
		SignalPeptideSensor.H \
		SignalPeptideSensor.C
	$(CC) $(CFLAGS) -o $(OBJ)/SignalPeptideSensor.o -c \
		SignalPeptideSensor.C
#---------------------------------------------------------
$(OBJ)/Isochore.o:\
		Isochore.H \
		Isochore.C
	$(CC) $(CFLAGS) -o $(OBJ)/Isochore.o -c \
		Isochore.C
#---------------------------------------------------------
#---------------------------------------------------------
$(OBJ)/SignalComparator.o:\
		SignalComparator.H \
		SignalComparator.C
	$(CC) $(CFLAGS) -o $(OBJ)/SignalComparator.o -c \
		SignalComparator.C
#---------------------------------------------------------
$(OBJ)/ModelBuilder.o:\
		genezilla.H \
		ModelBuilder.H \
		ModelBuilder.C
	$(CC) $(CFLAGS) -o $(OBJ)/ModelBuilder.o -c \
		ModelBuilder.C
#---------------------------------------------------------
$(OBJ)/SignalSensor.o:\
		genezilla.H \
		SignalSensor.H \
		SignalSensor.C
	$(CC) $(CFLAGS) -o $(OBJ)/SignalSensor.o -c \
		SignalSensor.C
#---------------------------------------------------------
$(OBJ)/SignalQueue.o:\
		SignalQueue.H \
		genezilla.H \
		SignalQueue.C \
		Propagator.H \
		Signal.H
	$(CC) $(CFLAGS) -o $(OBJ)/SignalQueue.o -c \
		SignalQueue.C
#---------------------------------------------------------
$(OBJ)/NthOrderStringIterator.o:\
		NthOrderStringIterator.H \
		NthOrderStringIterator.C
	$(CC) $(CFLAGS) -o $(OBJ)/NthOrderStringIterator.o -c \
		NthOrderStringIterator.C
#---------------------------------------------------------
$(OBJ)/MarkovChain.o:\
		genezilla.H \
		MarkovChain.H \
		MarkovChain.C
	$(CC) $(CFLAGS) -o $(OBJ)/MarkovChain.o -c \
		MarkovChain.C
#---------------------------------------------------------
$(OBJ)/ThreePeriodicIMM.o:\
		genezilla.H \
		ThreePeriodicIMM.H \
		ThreePeriodicIMM.C
	$(CC) $(CFLAGS) -o $(OBJ)/ThreePeriodicIMM.o -c \
		ThreePeriodicIMM.C
#---------------------------------------------------------
$(OBJ)/IMM.o:\
		genezilla.H \
		IMM.H \
		IMM.C
	$(CC) $(CFLAGS) -o $(OBJ)/IMM.o -c \
		IMM.C
#---------------------------------------------------------
$(OBJ)/ThreePeriodicMarkovChain.o:\
		genezilla.H \
		ThreePeriodicMarkovChain.H \
		ThreePeriodicMarkovChain.C
	$(CC) $(CFLAGS) -o $(OBJ)/ThreePeriodicMarkovChain.o -c \
		ThreePeriodicMarkovChain.C
#---------------------------------------------------------
$(OBJ)/ContentSensor.o:\
		ContentSensor.H \
		ContentSensor.C
	$(CC) $(CFLAGS) -o $(OBJ)/ContentSensor.o -c \
		ContentSensor.C
#---------------------------------------------------------
$(OBJ)/GarbageCollector.o:\
		GarbageCollector.H \
		GarbageCollector.C
	$(CC) $(CFLAGS) -o $(OBJ)/GarbageCollector.o -c \
		GarbageCollector.C
#---------------------------------------------------------
$(OBJ)/WMM.o:\
		WMM.H \
		WMM.C
	$(CC) $(CFLAGS) -o $(OBJ)/WMM.o -c \
		WMM.C
#---------------------------------------------------------
$(OBJ)/ModelType.o:\
		ModelType.H \
		ModelType.C
	$(CC) $(CFLAGS) -o $(OBJ)/ModelType.o -c \
		ModelType.C
#---------------------------------------------------------
$(OBJ)/ScoreAnalyzer.o:\
		ScoreAnalyzer.H \
		ScoreAnalyzer.C
	$(CC) $(CFLAGS) -o $(OBJ)/ScoreAnalyzer.o -c \
		ScoreAnalyzer.C
#---------------------------------------------------------
$(OBJ)/IntronQueue.o:\
		genezilla.H \
		IntronQueue.H \
		IntronQueue.C
	$(CC) $(CFLAGS) -o $(OBJ)/IntronQueue.o -c \
		IntronQueue.C
#---------------------------------------------------------
$(OBJ)/NoncodingQueue.o:\
		genezilla.H \
		NoncodingQueue.H \
		NoncodingQueue.C
	$(CC) $(CFLAGS) -o $(OBJ)/NoncodingQueue.o -c \
		NoncodingQueue.C
#---------------------------------------------------------
$(OBJ)/DiscreteDistribution.o:\
		DiscreteDistribution.H \
		DiscreteDistribution.C
	$(CC) $(CFLAGS) -o $(OBJ)/DiscreteDistribution.o -c \
		DiscreteDistribution.C
#---------------------------------------------------------
$(OBJ)/EmpiricalDistribution.o:\
		EmpiricalDistribution.H \
		EmpiricalDistribution.C
	$(CC) $(CFLAGS) -o $(OBJ)/EmpiricalDistribution.o -c \
		EmpiricalDistribution.C
#---------------------------------------------------------
$(OBJ)/GeometricDistribution.o:\
		GeometricDistribution.H \
		GeometricDistribution.C
	$(CC) $(CFLAGS) -o $(OBJ)/GeometricDistribution.o -c \
		GeometricDistribution.C
#---------------------------------------------------------
$(OBJ)/TataCapModel.o:\
		TataCapModel.H \
		TataCapModel.C
	$(CC) $(CFLAGS) -o $(OBJ)/TataCapModel.o -c \
		TataCapModel.C
#---------------------------------------------------------
$(OBJ)/TataCapSignal.o:\
		TataCapSignal.H \
		TataCapSignal.C
	$(CC) $(CFLAGS) -o $(OBJ)/TataCapSignal.o -c \
		TataCapSignal.C
#---------------------------------------------------------
$(OBJ)/Signal.o:\
		Signal.H \
		Signal.C
	$(CC) $(CFLAGS) -o $(OBJ)/Signal.o -c \
		Signal.C
#---------------------------------------------------------
$(OBJ)/SignalTypeProperties.o:\
		SignalTypeProperties.H \
		SignalTypeProperties.C
	$(CC) $(CFLAGS) -o $(OBJ)/SignalTypeProperties.o -c \
		SignalTypeProperties.C
#---------------------------------------------------------
$(OBJ)/SignalType.o:\
		SignalType.H \
		SignalType.C
	$(CC) $(CFLAGS) -o $(OBJ)/SignalType.o -c \
		SignalType.C
#---------------------------------------------------------
$(OBJ)/ContentType.o:\
		ContentType.H \
		ContentType.C
	$(CC) $(CFLAGS) -o $(OBJ)/ContentType.o -c \
		ContentType.C
#---------------------------------------------------------
$(OBJ)/TrainingSequence.o:\
		TrainingSequence.H \
		TrainingSequence.C
	$(CC) $(CFLAGS) -o $(OBJ)/TrainingSequence.o -c \
		TrainingSequence.C
#---------------------------------------------------------
$(OBJ)/Transitions.o:\
		genezilla.H \
		Transitions.H \
		Transitions.C
	$(CC) $(CFLAGS) -o $(OBJ)/Transitions.o -c \
		Transitions.C
#---------------------------------------------------------
$(OBJ)/MarkovChainCompiler.o:\
		MarkovChainCompiler.H \
		MarkovChainCompiler.C
	$(CC) $(CFLAGS) -o $(OBJ)/MarkovChainCompiler.o -c \
		MarkovChainCompiler.C
#---------------------------------------------------------
$(OBJ)/FastMarkovChain.o:\
		genezilla.H \
		FastMarkovChain.H \
		FastMarkovChain.C
	$(CC) $(CFLAGS) -o $(OBJ)/FastMarkovChain.o -c \
		FastMarkovChain.C
#---------------------------------------------------------
$(OBJ)/Fast3PMC.o:\
		genezilla.H \
		Fast3PMC.H \
		Fast3PMC.C
	$(CC) $(CFLAGS) -o $(OBJ)/Fast3PMC.o -c \
		Fast3PMC.C
#---------------------------------------------------------
$(OBJ)/genezilla.o:\
		genezilla.H \
		genezilla.C \
		Propagator.H \
		SignalQueue.H \
		$(BOOM)/Histogram.H \
		Signal.H
	$(CC) $(CFLAGS) -o $(OBJ)/genezilla.o -c \
		genezilla.C
#---------------------------------------------------------
$(OBJ)/GZilla.o:\
		$(BOOM)/Histogram.H \
		GZilla.H \
		GZilla.C
	$(CC) $(CFLAGS) -o $(OBJ)/GZilla.o -c \
		GZilla.C
#---------------------------------------------------------
$(OBJ)/crf.o:\
		crf.H \
		crf.C \
		Propagator.H \
		SignalQueue.H \
		$(BOOM)/Histogram.H \
		Signal.H
	$(CC) $(CFLAGS) -o $(OBJ)/crf.o -c \
		crf.C
#---------------------------------------------------------
$(OBJ)/CRF.o:\
		$(BOOM)/Histogram.H \
		CRF.H \
		CRF.C
	$(CC) $(CFLAGS) -o $(OBJ)/CRF.o -c \
		CRF.C
#---------------------------------------------------------
$(OBJ)/TopologyLoader.o:\
		TopologyLoader.H \
		TopologyLoader.C
	$(CC) $(CFLAGS) -o $(OBJ)/TopologyLoader.o -c \
		TopologyLoader.C
#---------------------------------------------------------
$(OBJ)/Edge.o:\
		Edge.H \
		Edge.C
	$(CC) $(CFLAGS) -o $(OBJ)/Edge.o -c \
		Edge.C
#---------------------------------------------------------
$(OBJ)/EdgeFactory.o:\
		EdgeFactory.H \
		EdgeFactory.C
	$(CC) $(CFLAGS) -o $(OBJ)/EdgeFactory.o -c \
		EdgeFactory.C
#---------------------------------------------------------
$(OBJ)/WAM.o:\
		WAM.H \
		WAM.C
	$(CC) $(CFLAGS) -o $(OBJ)/WAM.o -c \
		WAM.C
#---------------------------------------------------------
$(OBJ)/BranchAcceptor.o:\
		BranchAcceptor.H \
		BranchAcceptor.C
	$(CC) $(CFLAGS) -o $(OBJ)/BranchAcceptor.o -c \
		BranchAcceptor.C
#---------------------------------------------------------
$(OBJ)/ParseGraph.o:\
		ParseGraph.H \
		ParseGraph.C \
		Edge.H \
		Signal.H
	$(CC) $(CFLAGS) -o $(OBJ)/ParseGraph.o -c \
		ParseGraph.C
#---------------------------------------------------------
$(OBJ)/WWAM.o:\
		WWAM.H \
		WWAM.C
	$(CC) $(CFLAGS) -o $(OBJ)/WWAM.o -c \
		WWAM.C
#---------------------------------------------------------
$(OBJ)/reweight-graph.o:\
		reweight-graph.C
	$(CC) $(CFLAGS) -o $(OBJ)/reweight-graph.o -c \
		reweight-graph.C
#---------------------------------------------------------
reweight-graph: \
		$(OBJ)/reweight-graph.o \
		$(OBJ)/IntronDepthProfile.o \
		$(OBJ)/LightGraph.o \
		$(OBJ)/LightEdge.o \
		$(OBJ)/LightVertex.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/RnaJunction.o
	$(CC) $(LDFLAGS) -o reweight-graph \
		$(OBJ)/reweight-graph.o \
		$(OBJ)/IntronDepthProfile.o \
		$(OBJ)/LightGraph.o \
		$(OBJ)/LightEdge.o \
		$(OBJ)/LightVertex.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/RnaJunction.o \
		$(LIBDIRS) $(LIBS)
#---------------------------------------------------------
$(OBJ)/n-best.o:\
		n-best.C
	$(CC) $(CFLAGS) -o $(OBJ)/n-best.o -c \
		n-best.C
#---------------------------------------------------------
n-best: \
		$(OBJ)/NMD.o \
		$(OBJ)/n-best.o \
		$(OBJ)/LightGraph.o \
		$(OBJ)/LightEdge.o \
		$(OBJ)/LightVertex.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/TrellisLink.o \
		$(OBJ)/ContentType.o
	$(CC) $(LDFLAGS) -o n-best \
		$(OBJ)/NMD.o \
		$(OBJ)/n-best.o \
		$(OBJ)/LightGraph.o \
		$(OBJ)/LightEdge.o \
		$(OBJ)/TrellisLink.o \
		$(OBJ)/LightVertex.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(LIBDIRS) $(LIBS)
#---------------------------------------------------------
$(OBJ)/get-bias-model.o:\
		get-bias-model.C
	$(CC) $(CFLAGS) -o $(OBJ)/get-bias-model.o -c \
		get-bias-model.C
#---------------------------------------------------------
get-bias-model: \
		$(OBJ)/get-bias-model.o \
		$(OBJ)/IntronDepthProfile.o \
		$(OBJ)/LightGraph.o \
		$(OBJ)/LightEdge.o \
		$(OBJ)/LightVertex.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/RnaJunction.o
	$(CC) $(LDFLAGS) -o get-bias-model \
		$(OBJ)/get-bias-model.o \
		$(OBJ)/IntronDepthProfile.o \
		$(OBJ)/LightGraph.o \
		$(OBJ)/LightEdge.o \
		$(OBJ)/LightVertex.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/RnaJunction.o \
		$(LIBDIRS) $(LIBS)
#---------------------------------------------------------
genezilla: \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/genezilla.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/GZilla.o
	$(CC) $(LDFLAGS) -o genezilla \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/genezilla.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/GZilla.o \
		$(LIBDIRS) $(LIBS)
#---------------------------------------------------------
crf: \
		$(OBJ)/SignalLabelingProfile.o \
		$(OBJ)/LabelMatrix.o \
		$(OBJ)/Labeling.o \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/crf.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/GZilla.o \
		$(OBJ)/CRF.o
	$(CC) $(LDFLAGS) -o crf \
		$(OBJ)/SignalLabelingProfile.o \
		$(OBJ)/LabelMatrix.o \
		$(OBJ)/Labeling.o \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/crf.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/GZilla.o \
		$(OBJ)/CRF.o \
		$(LIBDIRS) $(LIBS)
#---------------------------------------------------------
genezilla-abinitio: \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/genezilla.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/GZilla.o
	$(CC) $(LDFLAGS) -o genezilla-abinitio \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/genezilla.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/GZilla.o \
		$(LIBDIRS) $(LIBS)
#---------------------------------------------------------
$(OBJ)/train-signal-sensor.o:\
		train-signal-sensor.C
	$(CC) $(CFLAGS) -o $(OBJ)/train-signal-sensor.o -c \
		train-signal-sensor.C
#---------------------------------------------------------
train-signal-sensor: \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/train-signal-sensor.o
	$(CC) $(LDFLAGS) -o train-signal-sensor \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Signal.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/train-signal-sensor.o \
		$(LIBDIRS) $(LIBS)

#---------------------------------------------------------
$(OBJ)/BWM.o:\
		BWM.H \
		BWM.C
	$(CC) $(CFLAGS) -o $(OBJ)/BWM.o -c \
		BWM.C
#---------------------------------------------------------
$(OBJ)/train-BWM.o:\
		train-BWM.C \
		BWM.H \
		BWM.C
	$(CC) $(CFLAGS) -o $(OBJ)/train-BWM.o -c \
		train-BWM.C
#---------------------------------------------------------
train-BWM: \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/BWM.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/train-BWM.o
	$(CC) $(LDFLAGS) -o train-BWM \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/BWM.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Signal.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/train-BWM.o \
		$(LIBDIRS) $(LIBS)

#---------------------------------------------------------
$(OBJ)/train-content-sensor.o:\
		train-content-sensor.C
	$(CC) $(CFLAGS) -o $(OBJ)/train-content-sensor.o -c \
		train-content-sensor.C
#---------------------------------------------------------
train-content-sensor: \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/Edge.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/train-content-sensor.o
	$(CC) $(LDFLAGS) -o train-content-sensor \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Edge.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Signal.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/train-content-sensor.o \
		$(LIBDIRS) $(LIBS)


#--------------------------------------------------------
$(OBJ)/compile-markov-chain.o:\
		compile-markov-chain.C
	$(CC) $(CFLAGS) -o $(OBJ)/compile-markov-chain.o -c \
		compile-markov-chain.C
#---------------------------------------------------------
compile-markov-chain: \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/compile-markov-chain.o
	$(CC) $(LDFLAGS) -o compile-markov-chain \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/Edge.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/compile-markov-chain.o \
		$(LIBDIRS) $(LIBS)

#---------------------------------------------
$(OBJ)/mdd.o:\
		mdd.C
	$(CC) $(CFLAGS) -o $(OBJ)/mdd.o -c \
		mdd.C
#---------------------------------------------------------
$(OBJ)/MddTree.o:\
		MddTree.C \
		MddTree.H
	$(CC) $(CFLAGS) -o $(OBJ)/MddTree.o -c \
		MddTree.C
#---------------------------------------------------------
$(OBJ)/Partition.o:\
		Partition.C \
		Partition.H
	$(CC) $(CFLAGS) -o $(OBJ)/Partition.o -c \
		Partition.C
#---------------------------------------------------------
$(OBJ)/TreeNode.o:\
		TreeNode.C \
		TreeNode.H
	$(CC) $(CFLAGS) -o $(OBJ)/TreeNode.o -c \
		TreeNode.C
#---------------------------------------------------------
mdd: \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/Partition.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/mdd.o
	$(CC) $(LDFLAGS) -o mdd \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Signal.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/mdd.o \
		$(LIBDIRS) $(LIBS)

#---------------------------------------------------------



#--------------------------------------------------------
$(OBJ)/refine-polya.o:\
		refine-polya.C
	$(CC) $(CFLAGS) -o $(OBJ)/refine-polya.o -c \
		refine-polya.C
#---------------------------------------------------------
refine-polya: \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/refine-polya.o
	$(CC) $(LDFLAGS) -o refine-polya \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/refine-polya.o \
		$(LIBDIRS) $(LIBS)

#---------------------------------------------



#--------------------------------------------------------
$(OBJ)/refine-tata.o:\
		refine-tata.C
	$(CC) $(CFLAGS) -o $(OBJ)/refine-tata.o -c \
		refine-tata.C
#---------------------------------------------------------
refine-tata: \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/refine-tata.o
	$(CC) $(LDFLAGS) -o refine-tata \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/refine-tata.o \
		$(LIBDIRS) $(LIBS)

#---------------------------------------------------------
$(OBJ)/get-classifier-examples.o:\
		get-classifier-examples.C \
		GZilla.C \
		GZilla.H
	$(CC) $(CFLAGS) -o $(OBJ)/get-classifier-examples.o -c \
		get-classifier-examples.C
#---------------------------------------------------------
get-classifier-examples: \
		$(OBJ)/get-classifier-examples.o\
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Signal.o \
		$(OBJ)/GZilla.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/SignalTypeProperties.o \
	$(CC) $(LDFLAGS) -o get-classifier-examples \
		$(OBJ)/get-classifier-examples.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/GZilla.o \
		$(OBJ)/Signal.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/SignalTypeProperties.o \
		$(LIBDIRS) $(LIBS)

#--------------------------------------------------------
$(OBJ)/orf-stats.o:\
		orf-stats.C
	$(CC) $(CFLAGS) -o $(OBJ)/orf-stats.o -c \
		orf-stats.C
#---------------------------------------------------------
orf-stats: \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Signal.o \
		$(OBJ)/GZilla.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/orf-stats.o
	$(CC) $(LDFLAGS) -o orf-stats \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Signal.o \
		$(OBJ)/GZilla.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/orf-stats.o \
		$(STATIC) \
		$(LIBDIRS) $(LIBS)
#--------------------------------------------------------
$(OBJ)/entropy.o:\
		entropy.C
	$(CC) $(CFLAGS) -o $(OBJ)/entropy.o -c \
		entropy.C
#---------------------------------------------------------
entropy: \
		$(OBJ)/entropy.o
	$(CC) $(LDFLAGS) -o entropy \
		$(OBJ)/entropy.o \
		$(STATIC)
#---------------------------------------------
#---------------------------------------------------------
$(OBJ)/score-gff.o:\
		genezilla.H \
		score-gff.C \
		Propagator.H \
		SignalQueue.H \
		Signal.H
	$(CC) $(CFLAGS) -o $(OBJ)/score-gff.o -c \
		score-gff.C
#---------------------------------------------------------
score-gff: \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/score-gff.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/GZilla.o \
	$(CC) $(LDFLAGS) -o score-gff \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/score-gff.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/GZilla.o \
		$(STATIC) \
		$(LIBDIRS) $(LIBS)

#--------------------------------------------------------
$(OBJ)/WMM-add-pseudocounts.o:\
		WMM-add-pseudocounts.C
	$(CC) $(CFLAGS) -o $(OBJ)/WMM-add-pseudocounts.o -c \
		WMM-add-pseudocounts.C
#---------------------------------------------------------
WMM-add-pseudocounts: \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/Partition.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/WMM-add-pseudocounts.o
	$(CC) $(LDFLAGS) -o WMM-add-pseudocounts \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/Partition.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/WMM-add-pseudocounts.o \
		$(STATIC) \
		$(LIBDIRS) $(LIBS)

#---------------------------------------------

#--------------------------------------------------------
$(OBJ)/compute-signal-likelihoods.o:\
		compute-signal-likelihoods.C \
		$(BOOM)/Histogram.H
	$(CC) $(CFLAGS) -o $(OBJ)/compute-signal-likelihoods.o -c \
		compute-signal-likelihoods.C
#---------------------------------------------------------
compute-signal-likelihoods: \
		$(OBJ)/compute-signal-likelihoods.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/GZilla.o
	$(CC) $(LDFLAGS) -o compute-signal-likelihoods \
		$(OBJ)/compute-signal-likelihoods.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/GZilla.o \
		$(LIBS)

#---------------------------------------------

#--------------------------------------------------------
$(OBJ)/train-branch-acceptor.o:\
		train-branch-acceptor.C
	$(CC) $(CFLAGS) -o $(OBJ)/train-branch-acceptor.o -c \
		train-branch-acceptor.C
#---------------------------------------------------------
train-branch-acceptor: \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Signal.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/train-branch-acceptor.o
	$(CC) $(LDFLAGS) -o train-branch-acceptor \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Signal.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/train-branch-acceptor.o \
		$(STATIC) \
		$(LIBDIRS) $(LIBS)

#---------------------------------------------

#--------------------------------------------------------
$(OBJ)/train-signal-peptide-model.o:\
		train-signal-peptide-model.C
	$(CC) $(CFLAGS) -o $(OBJ)/train-signal-peptide-model.o -c \
		train-signal-peptide-model.C
#---------------------------------------------------------
train-signal-peptide-model: \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/IMM.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/Partition.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/WMM.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Edge.o \
		$(OBJ)/train-signal-peptide-model.o
	$(CC) $(LDFLAGS) -o train-signal-peptide-model \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/Edge.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/IMM.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/Partition.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/WMM.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/train-signal-peptide-model.o \
		$(LIBDIRS) $(LIBS)

#---------------------------------------------

#--------------------------------------------------------
$(OBJ)/train-tata-cap-model.o:\
		train-tata-cap-model.C \
		train-tata-cap-model.H
	$(CC) $(CFLAGS) -o $(OBJ)/train-tata-cap-model.o -c \
		train-tata-cap-model.C
#---------------------------------------------------------
train-tata-cap-model: \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Signal.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/train-tata-cap-model.o
	$(CC) $(LDFLAGS) -o train-tata-cap-model \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/IMM.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Signal.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/train-tata-cap-model.o \
		$(LIBDIRS) $(LIBS)

#---------------------------------------------

#--------------------------------------------------------
$(OBJ)/find-cpg-islands.o:\
		find-cpg-islands.C
	$(CC) $(CFLAGS) -o $(OBJ)/find-cpg-islands.o -c \
		find-cpg-islands.C
#---------------------------------------------------------
find-cpg-islands: \
		$(OBJ)/find-cpg-islands.o
	$(CC) $(LDFLAGS) -o find-cpg-islands \
		$(OBJ)/find-cpg-islands.o
#---------------------------------------------

#--------------------------------------------------------
$(OBJ)/get-transcripts.o:\
		get-transcripts.C
	$(CC) $(CFLAGS) -o $(OBJ)/get-transcripts.o -c \
		get-transcripts.C
#---------------------------------------------------------
get-transcripts: \
		$(OBJ)/get-transcripts.o
	$(CC) $(LDFLAGS) -o get-transcripts \
		$(OBJ)/get-transcripts.o \
		$(LIBDIRS) $(LIBS)
#---------------------------------------------


#---------------------------------------------------------
$(OBJ)/apply-content-sensor.o:\
		apply-content-sensor.C
	$(CC) $(CFLAGS) -o $(OBJ)/apply-content-sensor.o -c \
		apply-content-sensor.C
#---------------------------------------------------------
apply-content-sensor: \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/Edge.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/apply-content-sensor.o
	$(CC) $(LDFLAGS) -o apply-content-sensor \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Edge.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Signal.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/apply-content-sensor.o \
		$(LIBDIRS) $(LIBS)
#---------------------------------------------------------
$(OBJ)/apply-signal-sensor.o:\
		apply-signal-sensor.C
	$(CC) $(CFLAGS) -o $(OBJ)/apply-signal-sensor.o -c \
		apply-signal-sensor.C
#---------------------------------------------------------
apply-signal-sensor: \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/Edge.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/apply-signal-sensor.o
	$(CC) $(LDFLAGS) -o apply-signal-sensor \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Edge.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/Signal.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/apply-signal-sensor.o \
		$(LIBDIRS) $(LIBS)
#---------------------------------------------------------
$(OBJ)/fit.o:\
		fit.C
	$(CC) $(CFLAGS) -o $(OBJ)/fit.o -c \
		fit.C
#---------------------------------------------------------
fit: \
		$(OBJ)/fit.o
	$(CC) $(LDFLAGS) -o fit \
		$(OBJ)/fit.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/score-content.o:\
		score-content.C
	$(CC) $(CFLAGS) -o $(OBJ)/score-content.o -c \
		score-content.C
#---------------------------------------------------------
score-content: \
		$(OBJ)/score-content.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/Edge.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o
	$(CC) $(LDFLAGS) -o score-content \
		$(OBJ)/score-content.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/Edge.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/relative-entropy.o:\
		relative-entropy.C
	$(CC) $(CFLAGS) -o $(OBJ)/relative-entropy.o -c \
		relative-entropy.C
#---------------------------------------------------------
relative-entropy: \
		$(OBJ)/relative-entropy.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/Edge.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o
	$(CC) $(LDFLAGS) -o relative-entropy \
		$(OBJ)/relative-entropy.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/Edge.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/score-sequence.o:\
		score-sequence.C
	$(CC) $(CFLAGS) -o $(OBJ)/score-sequence.o -c \
		score-sequence.C
#---------------------------------------------------------
score-sequence: \
		$(OBJ)/score-sequence.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/Edge.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o
	$(CC) $(LDFLAGS) -o score-sequence \
		$(OBJ)/score-sequence.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/Edge.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/classify-content.o:\
		classify-content.C
	$(CC) $(CFLAGS) -o $(OBJ)/classify-content.o -c \
		classify-content.C
#---------------------------------------------------------
classify-content: \
		$(OBJ)/classify-content.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/Edge.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o
	$(CC) $(LDFLAGS) -o classify-content \
		$(OBJ)/classify-content.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/Edge.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/WMM.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/ModelType.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(LIBS)
#--------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/EvidenceFilter.o:\
		EvidenceFilter.C\
		EvidenceFilter.H
	$(CC) $(CFLAGS) -o $(OBJ)/EvidenceFilter.o -c \
		EvidenceFilter.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/LightGraph.o:\
		LightGraph.C\
		LightGraph.H
	$(CC) $(CFLAGS) -o $(OBJ)/LightGraph.o -c \
		LightGraph.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/LightEdge.o:\
		LightEdge.C\
		LightEdge.H
	$(CC) $(CFLAGS) -o $(OBJ)/LightEdge.o -c \
		LightEdge.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/LightVertex.o:\
		LightVertex.C\
		LightVertex.H
	$(CC) $(CFLAGS) -o $(OBJ)/LightVertex.o -c \
		LightVertex.C
#---------------------------------------------------------
#--------------------------------------------------------
$(OBJ)/IntronDepthProfile.o:\
		IntronDepthProfile.C\
		IntronDepthProfile.H
	$(CC) $(CFLAGS) -o $(OBJ)/IntronDepthProfile.o -c \
		IntronDepthProfile.C
#---------------------------------------------------------
$(OBJ)/TrellisLink.o:\
		TrellisLink.C\
		TrellisLink.H
	$(CC) $(CFLAGS) -o $(OBJ)/TrellisLink.o -c \
		TrellisLink.C
#--------------------------------------------------------
$(OBJ)/project-annotation.o:\
		project-annotation.C
	$(CC) $(CFLAGS) -o $(OBJ)/project-annotation.o -c \
		project-annotation.C
#---------------------------------------------------------
project-annotation: \
		$(OBJ)/Labeling.o \
		$(OBJ)/ProjectionChecker.o \
		$(OBJ)/project-annotation.o
	$(CC) $(LDFLAGS) -o project-annotation \
		$(OBJ)/Labeling.o \
		$(OBJ)/ProjectionChecker.o \
		$(OBJ)/project-annotation.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/Labeling.o:\
		Labeling.C\
		Labeling.H
	$(CC) $(CFLAGS) -o $(OBJ)/Labeling.o -c \
		Labeling.C
#--------------------------------------------------------
$(OBJ)/LabelMatrix.o:\
		LabelMatrix.C\
		LabelMatrix.H
	$(CC) $(CFLAGS) -o $(OBJ)/LabelMatrix.o -c \
		LabelMatrix.C
#--------------------------------------------------------
$(OBJ)/SignalLabelingProfile.o:\
		SignalLabelingProfile.C\
		SignalLabelingProfile.H
	$(CC) $(CFLAGS) -o $(OBJ)/SignalLabelingProfile.o -c \
		SignalLabelingProfile.C
#--------------------------------------------------------
$(OBJ)/posmod.o:\
		posmod.C
	$(CC) $(CFLAGS) -o $(OBJ)/posmod.o -c \
		posmod.C
#---------------------------------------------------------
posmod: \
		$(OBJ)/posmod.o
	$(CC) $(LDFLAGS) -o posmod \
		$(OBJ)/posmod.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/AlignmentPath.o:\
		AlignmentPath.H \
		AlignmentPath.C
	$(CC) $(CFLAGS) -o $(OBJ)/AlignmentPath.o -c \
		AlignmentPath.C
#---------------------------------------------------------
$(OBJ)/align.o:\
		align.C
	$(CC) $(CFLAGS) -o $(OBJ)/align.o -c \
		align.C
#---------------------------------------------------------
align: \
		$(OBJ)/AlignmentPath.o \
		$(OBJ)/align.o
	$(CC) $(LDFLAGS) -o align \
		$(OBJ)/AlignmentPath.o \
		$(OBJ)/align.o \
		$(LIBS)
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/find-variants.o:\
		find-variants.C
	$(CC) $(CFLAGS) -o $(OBJ)/find-variants.o -c \
		find-variants.C
#---------------------------------------------------------
find-variants: \
		$(OBJ)/find-variants.o
	$(CC) $(LDFLAGS) -o find-variants \
		$(OBJ)/find-variants.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/pds.o:\
		pds.C
	$(CC) $(CFLAGS) -o $(OBJ)/pds.o -c \
		pds.C
#---------------------------------------------------------
pds: \
		$(OBJ)/pds.o
	$(CC) $(LDFLAGS) -o pds \
		$(OBJ)/pds.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/check-projection.o:\
		check-projection.C
	$(CC) $(CFLAGS) -o $(OBJ)/check-projection.o -c \
		check-projection.C
#---------------------------------------------------------
check-projection: \
		$(OBJ)/Labeling.o \
		$(OBJ)/ProjectionChecker.o \
		$(OBJ)/check-projection.o
	$(CC) $(LDFLAGS) -o check-projection \
		$(OBJ)/Labeling.o \
		$(OBJ)/ProjectionChecker.o \
		$(OBJ)/check-projection.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/vcf-to-gcf.o:\
		vcf-to-gcf.C
	$(CC) $(CFLAGS) -o $(OBJ)/vcf-to-gcf.o -c \
		vcf-to-gcf.C
#---------------------------------------------------------
vcf-to-gcf: \
		$(OBJ)/vcf-to-gcf.o
	$(CC) $(LDFLAGS) -o vcf-to-gcf \
		$(OBJ)/vcf-to-gcf.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/gcf-to-fasta.o:\
		gcf-to-fasta.C
	$(CC) $(CFLAGS) -o $(OBJ)/gcf-to-fasta.o -c \
		gcf-to-fasta.C
#---------------------------------------------------------
gcf-to-fasta: \
		$(OBJ)/gcf-to-fasta.o
	$(CC) $(LDFLAGS) -o gcf-to-fasta \
		$(OBJ)/gcf-to-fasta.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/gcf-list-samples.o:\
		gcf-list-samples.C
	$(CC) $(CFLAGS) -o $(OBJ)/gcf-list-samples.o -c \
		gcf-list-samples.C
#---------------------------------------------------------
gcf-list-samples: \
		$(OBJ)/gcf-list-samples.o
	$(CC) $(LDFLAGS) -o gcf-list-samples \
		$(OBJ)/gcf-list-samples.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/split-fasta-by-tag.o:\
		split-fasta-by-tag.C
	$(CC) $(CFLAGS) -o $(OBJ)/split-fasta-by-tag.o -c \
		split-fasta-by-tag.C
#---------------------------------------------------------
split-fasta-by-tag: \
		$(OBJ)/split-fasta-by-tag.o
	$(CC) $(LDFLAGS) -o split-fasta-by-tag \
		$(OBJ)/split-fasta-by-tag.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/split-gff-by-gene.o:\
		split-gff-by-gene.C
	$(CC) $(CFLAGS) -o $(OBJ)/split-gff-by-gene.o -c \
		split-gff-by-gene.C
#---------------------------------------------------------
split-gff-by-gene: \
		$(OBJ)/split-gff-by-gene.o
	$(CC) $(LDFLAGS) -o split-gff-by-gene \
		$(OBJ)/split-gff-by-gene.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/gff-fix-chrom-names.o:\
		gff-fix-chrom-names.C
	$(CC) $(CFLAGS) -o $(OBJ)/gff-fix-chrom-names.o -c \
		gff-fix-chrom-names.C
#---------------------------------------------------------
gff-fix-chrom-names: \
		$(OBJ)/gff-fix-chrom-names.o
	$(CC) $(LDFLAGS) -o gff-fix-chrom-names \
		$(OBJ)/gff-fix-chrom-names.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/gff-shift-coords.o:\
		gff-shift-coords.C
	$(CC) $(CFLAGS) -o $(OBJ)/gff-shift-coords.o -c \
		gff-shift-coords.C
#---------------------------------------------------------
gff-shift-coords: \
		$(OBJ)/gff-shift-coords.o
	$(CC) $(LDFLAGS) -o gff-shift-coords \
		$(OBJ)/gff-shift-coords.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/fasta-get-id.o:\
		fasta-get-id.C
	$(CC) $(CFLAGS) -o $(OBJ)/fasta-get-id.o -c \
		fasta-get-id.C
#---------------------------------------------------------
fasta-get-id: \
		$(OBJ)/fasta-get-id.o
	$(CC) $(LDFLAGS) -o fasta-get-id \
		$(OBJ)/fasta-get-id.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/gff-change-substrate.o:\
		gff-change-substrate.C
	$(CC) $(CFLAGS) -o $(OBJ)/gff-change-substrate.o -c \
		gff-change-substrate.C
#---------------------------------------------------------
gff-change-substrate: \
		$(OBJ)/gff-change-substrate.o
	$(CC) $(LDFLAGS) -o gff-change-substrate \
		$(OBJ)/gff-change-substrate.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/fastas-are-identical.o:\
		fastas-are-identical.C
	$(CC) $(CFLAGS) -o $(OBJ)/fastas-are-identical.o -c \
		fastas-are-identical.C
#---------------------------------------------------------
fastas-are-identical: \
		$(OBJ)/fastas-are-identical.o
	$(CC) $(LDFLAGS) -o fastas-are-identical \
		$(OBJ)/fastas-are-identical.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/gene-LD.o:\
		gene-LD.C
	$(CC) $(CFLAGS) -o $(OBJ)/gene-LD.o -c \
		gene-LD.C
#---------------------------------------------------------
gene-LD: \
		$(OBJ)/gene-LD.o
	$(CC) $(LDFLAGS) -o gene-LD \
		$(OBJ)/gene-LD.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/gff-longest-transcript.o:\
		gff-longest-transcript.C
	$(CC) $(CFLAGS) -o $(OBJ)/gff-longest-transcript.o -c \
		gff-longest-transcript.C
#---------------------------------------------------------
gff-longest-transcript: \
		$(OBJ)/gff-longest-transcript.o
	$(CC) $(LDFLAGS) -o gff-longest-transcript \
		$(OBJ)/gff-longest-transcript.o \
		$(LIBS)
#---------------------------------------------------------
find-variant-signals: \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/find-variant-signals.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/GZilla.o
	$(CC) $(LDFLAGS) -o find-variant-signals \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/find-variant-signals.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/GZilla.o \
		$(LIBDIRS) $(LIBS)
#--------------------------------------------------------
$(OBJ)/find-variant-signals.o:\
		find-variant-signals.C
	$(CC) $(CFLAGS) -o $(OBJ)/find-variant-signals.o -c \
		find-variant-signals.C
#---------------------------------------------------------
$(OBJ)/ReferenceAnnotation.o:\
		ReferenceAnnotation.C\
		ReferenceAnnotation.H
	$(CC) $(CFLAGS) -o $(OBJ)/ReferenceAnnotation.o -c \
		ReferenceAnnotation.C
#---------------------------------------------------------
$(OBJ)/ContentRegion.o:\
		ContentRegion.C\
		ContentRegion.H
	$(CC) $(CFLAGS) -o $(OBJ)/ContentRegion.o -c \
		ContentRegion.C
#---------------------------------------------------------
$(OBJ)/ContentRegions.o:\
		ContentRegions.C\
		ContentRegions.H
	$(CC) $(CFLAGS) -o $(OBJ)/ContentRegions.o -c \
		ContentRegions.C
#---------------------------------------------------------
$(OBJ)/VariantEvent.o:\
		VariantEvent.C\
		VariantEvent.H
	$(CC) $(CFLAGS) -o $(OBJ)/VariantEvent.o -c \
		VariantEvent.C
#---------------------------------------------------------
$(OBJ)/VariantEvents.o:\
		VariantEvents.C\
		VariantEvents.H
	$(CC) $(CFLAGS) -o $(OBJ)/VariantEvents.o -c \
		VariantEvents.C
#---------------------------------------------------------
$(OBJ)/SignalStream.o:\
		SignalStream.C\
		SignalStream.H
	$(CC) $(CFLAGS) -o $(OBJ)/SignalStream.o -c \
		SignalStream.C
#---------------------------------------------------------
$(OBJ)/ConstraintInterval.o:\
		ConstraintInterval.C\
		ConstraintInterval.H
	$(CC) $(CFLAGS) -o $(OBJ)/ConstraintInterval.o -c \
		ConstraintInterval.C
#---------------------------------------------------------
$(OBJ)/ConstraintIntervals.o:\
		ConstraintIntervals.C\
		ConstraintIntervals.H
	$(CC) $(CFLAGS) -o $(OBJ)/ConstraintIntervals.o -c \
		ConstraintIntervals.C
#---------------------------------------------------------
$(OBJ)/SignalStreamBuilder.o:\
		SignalStreamBuilder.C \
		SignalStreamBuilder.H \
		ContentRegion.H \
		ContentRegions.H \
		ReferenceAnnotation.H \
		VariantEvent.H \
		VariantEvents.H \
		SignalStream.H \
		ConstraintInterval.H \
		ConstraintIntervals.H
	$(CC) $(CFLAGS) -o $(OBJ)/SignalStreamBuilder.o -c \
		SignalStreamBuilder.C
#---------------------------------------------------------
$(OBJ)/cia.o:\
		SignalStreamBuilder.H \
		ContentRegion.H \
		ContentRegions.H \
		ReferenceAnnotation.H \
		VariantEvent.H \
		VariantEvents.H \
		SignalStream.H \
		ConstraintInterval.H \
		ConstraintIntervals.H \
		cia.H \
		cia.C \
		Propagator.H \
		SignalQueue.H \
		$(BOOM)/Histogram.H \
		Signal.H
	$(CC) $(CFLAGS) -o $(OBJ)/cia.o -c \
		cia.C
#---------------------------------------------------------
$(OBJ)/CIA.o:\
		$(BOOM)/Histogram.H \
		CIA.H \
		CIA.C
	$(CC) $(CFLAGS) -o $(OBJ)/CIA.o -c \
		CIA.C
#---------------------------------------------------------
cia: \
		$(OBJ)/SignalStreamBuilder.o \
		$(OBJ)/ContentRegion.o \
		$(OBJ)/ContentRegions.o \
		$(OBJ)/ReferenceAnnotation.o \
		$(OBJ)/VariantEvent.o \
		$(OBJ)/VariantEvents.o \
		$(OBJ)/SignalStream.o \
		$(OBJ)/ConstraintInterval.o \
		$(OBJ)/ConstraintIntervals.o \
		$(OBJ)/SignalLabelingProfile.o \
		$(OBJ)/LabelMatrix.o \
		$(OBJ)/Labeling.o \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/cia.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/GZilla.o \
		$(OBJ)/PriorMask.o \
		$(OBJ)/CIA.o
	$(CC) $(LDFLAGS) -o cia \
		$(OBJ)/PriorMask.o \
		$(OBJ)/SignalStreamBuilder.o \
		$(OBJ)/ContentRegion.o \
		$(OBJ)/ContentRegions.o \
		$(OBJ)/ReferenceAnnotation.o \
		$(OBJ)/VariantEvent.o \
		$(OBJ)/VariantEvents.o \
		$(OBJ)/SignalStream.o \
		$(OBJ)/ConstraintInterval.o \
		$(OBJ)/ConstraintIntervals.o \
		$(OBJ)/SignalLabelingProfile.o \
		$(OBJ)/LabelMatrix.o \
		$(OBJ)/Labeling.o \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/cia.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/GZilla.o \
		$(OBJ)/CIA.o \
		$(LIBDIRS) $(LIBS)
#---------------------------------------------------------
$(OBJ)/SpliceFeatureExtractor.o:\
		SpliceFeatureExtractor.C\
		SpliceFeatureExtractor.H
	$(CC) $(CFLAGS) -o $(OBJ)/SpliceFeatureExtractor.o -c \
		SpliceFeatureExtractor.C
#---------------------------------------------------------
$(OBJ)/extract-splice-features.o:\
		extract-splice-features.C
	$(CC) $(CFLAGS) -o $(OBJ)/extract-splice-features.o -c \
		extract-splice-features.C
#---------------------------------------------------------
extract-splice-features: \
		$(OBJ)/SpliceFeatureExtractor.o \
		$(OBJ)/extract-splice-features.o \
		$(OBJ)/SignalStreamBuilder.o \
		$(OBJ)/ContentRegion.o \
		$(OBJ)/ContentRegions.o \
		$(OBJ)/ReferenceAnnotation.o \
		$(OBJ)/VariantEvent.o \
		$(OBJ)/VariantEvents.o \
		$(OBJ)/SignalStream.o \
		$(OBJ)/ConstraintInterval.o \
		$(OBJ)/ConstraintIntervals.o \
		$(OBJ)/SignalLabelingProfile.o \
		$(OBJ)/LabelMatrix.o \
		$(OBJ)/Labeling.o \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/GZilla.o \
		$(OBJ)/CIA.o
	$(CC) $(LDFLAGS) -o extract-splice-features \
		$(OBJ)/SpliceFeatureExtractor.o \
		$(OBJ)/SignalStreamBuilder.o \
		$(OBJ)/ContentRegion.o \
		$(OBJ)/ContentRegions.o \
		$(OBJ)/ReferenceAnnotation.o \
		$(OBJ)/VariantEvent.o \
		$(OBJ)/VariantEvents.o \
		$(OBJ)/SignalStream.o \
		$(OBJ)/ConstraintInterval.o \
		$(OBJ)/ConstraintIntervals.o \
		$(OBJ)/SignalLabelingProfile.o \
		$(OBJ)/LabelMatrix.o \
		$(OBJ)/Labeling.o \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/GZilla.o \
		$(OBJ)/CIA.o \
		$(OBJ)/extract-splice-features.o \
		$(LIBDIRS) $(LIBS)
#--------------------------------------------------------
$(OBJ)/simulate-cryptic-sites.o:\
		simulate-cryptic-sites.C
	$(CC) $(CFLAGS) -o $(OBJ)/simulate-cryptic-sites.o -c \
		simulate-cryptic-sites.C
#---------------------------------------------------------
simulate-cryptic-sites: \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/Labeling.o \
		$(OBJ)/GCcontent.o \
		$(OBJ)/GZilla.o \
		$(OBJ)/simulate-cryptic-sites.o
	$(CC) $(LDFLAGS) -o simulate-cryptic-sites \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/Labeling.o \
		$(OBJ)/GCcontent.o \
		$(OBJ)/GZilla.o \
		$(OBJ)/simulate-cryptic-sites.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/GCcontent.o:\
		GCcontent.C\
		GCcontent.H
	$(CC) $(CFLAGS) -o $(OBJ)/GCcontent.o -c \
		GCcontent.C
#---------------------------------------------------------
$(OBJ)/PriorMask.o:\
		PriorMask.C\
		PriorMask.H
	$(CC) $(CFLAGS) -o $(OBJ)/PriorMask.o -c \
		PriorMask.C
#---------------------------------------------------------
$(OBJ)/NMD.o:\
		NMD.C\
		NMD.H
	$(CC) $(CFLAGS) -o $(OBJ)/NMD.o -c \
		NMD.C
#---------------------------------------------------------
$(OBJ)/ProjectionChecker.o:\
		ProjectionChecker.C\
		ProjectionChecker.H
	$(CC) $(CFLAGS) -o $(OBJ)/ProjectionChecker.o -c \
		ProjectionChecker.C
#---------------------------------------------------------
$(OBJ)/find-non-splice-sites.o:\
		find-non-splice-sites.C
	$(CC) $(CFLAGS) -o $(OBJ)/find-non-splice-sites.o -c \
		find-non-splice-sites.C
#---------------------------------------------------------
find-non-splice-sites: \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/Labeling.o \
		$(OBJ)/GCcontent.o \
		$(OBJ)/GZilla.o \
		$(OBJ)/find-non-splice-sites.o
	$(CC) $(LDFLAGS) -o find-non-splice-sites \
		$(OBJ)/EvidenceFilter.o \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/ParseGraph.o \
		$(OBJ)/GffPathFromParseGraph.o \
		$(OBJ)/SignalComparator.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/TrainingSequence.o \
		$(OBJ)/SignalPeptideSensor.o \
		$(OBJ)/CodonTree.o \
		$(OBJ)/Isochore.o \
		$(OBJ)/IsochoreTable.o \
		$(OBJ)/BranchAcceptor.o \
		$(OBJ)/ThreePeriodicIMM.o \
		$(OBJ)/IMM.o \
		$(OBJ)/EdgeFactory.o \
		$(OBJ)/MddTree.o \
		$(OBJ)/Partition.o \
		$(OBJ)/TreeNode.o \
		$(OBJ)/GarbageCollector.o \
		$(OBJ)/Edge.o \
		$(OBJ)/TopologyLoader.o \
		$(OBJ)/WAM.o \
		$(OBJ)/WWAM.o \
		$(OBJ)/TataCapModel.o \
		$(OBJ)/TataCapSignal.o \
		$(OBJ)/MarkovChainCompiler.o \
		$(OBJ)/Fast3PMC.o \
		$(OBJ)/FastMarkovChain.o \
		$(OBJ)/ThreePeriodicMarkovChain.o \
		$(OBJ)/DiscreteDistribution.o \
		$(OBJ)/Transitions.o \
		$(OBJ)/EmpiricalDistribution.o \
		$(OBJ)/GeometricDistribution.o \
		$(OBJ)/NoncodingQueue.o \
		$(OBJ)/IntronQueue.o \
		$(OBJ)/SignalType.o \
		$(OBJ)/ContentType.o \
		$(OBJ)/ModelBuilder.o \
		$(OBJ)/ScoreAnalyzer.o \
		$(OBJ)/ContentSensor.o \
		$(OBJ)/MarkovChain.o \
		$(OBJ)/WMM.o \
		$(OBJ)/SignalQueue.o \
		$(OBJ)/SignalSensor.o \
		$(OBJ)/Propagator.o \
		$(OBJ)/Signal.o \
		$(OBJ)/SignalTypeProperties.o \
		$(OBJ)/Labeling.o \
		$(OBJ)/GCcontent.o \
		$(OBJ)/GZilla.o \
		$(OBJ)/find-non-splice-sites.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/split-vcf-by-pop.o:\
		split-vcf-by-pop.C
	$(CC) $(CFLAGS) -o $(OBJ)/split-vcf-by-pop.o -c \
		split-vcf-by-pop.C
#---------------------------------------------------------
split-vcf-by-pop: \
		$(OBJ)/split-vcf-by-pop.o
	$(CC) $(LDFLAGS) -o split-vcf-by-pop \
		$(OBJ)/split-vcf-by-pop.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/remove-unpaired-reads.o:\
		remove-unpaired-reads.C
	$(CC) $(CFLAGS) -o $(OBJ)/remove-unpaired-reads.o -c \
		remove-unpaired-reads.C
#---------------------------------------------------------
remove-unpaired-reads: \
		$(OBJ)/remove-unpaired-reads.o
	$(CC) $(LDFLAGS) -o remove-unpaired-reads \
		$(OBJ)/remove-unpaired-reads.o \
		$(LIBS)
#----------------------------------------------------------
