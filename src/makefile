# LIBRARY SETTINGS - SET AS NECESSARY

LIB_BOOST = /Net/fs1/home/dilthey/PnP/libs/boost_1_41_0/
INCS = -I$(LIB_BOOST) -IGraph 

# END LIBRARY SETTINGS

#
# object and binary dirs
#

DIR_OBJ = ../obj
DIR_BIN = ../bin

CXX    = g++
COPTS  = -ggdb -O2 -fopenmp -std=gnu++0x -fstack-protector-all
CFLAGS = 
COMPILE = $(CXX) $(INCS) $(CFLAGS) $(COPTS)

VPATH = Data:Graph:NextGen:hash:hash/deBruijn:hash/sequence:GraphAligner:GraphAlignerUnique
        
OBJS = \
        $(DIR_OBJ)/GenotypePanel.o \
        $(DIR_OBJ)/HaplotypePanel.o \
        $(DIR_OBJ)/Edge.o \
        $(DIR_OBJ)/Graph.o \
        $(DIR_OBJ)/LargeGraph.o \
        $(DIR_OBJ)/MultiGraph.o \
        $(DIR_OBJ)/HMM.o \
        $(DIR_OBJ)/MultiHMM.o \
        $(DIR_OBJ)/AlphaHMM.o \
        $(DIR_OBJ)/AlphaHMMlegacy.o \
        $(DIR_OBJ)/Node.o \
        $(DIR_OBJ)/NextGen.o \
        $(DIR_OBJ)/simulationSuite.o \
        $(DIR_OBJ)/Validation.o \
        $(DIR_OBJ)/LocusCodeAllocation.o \
        $(DIR_OBJ)/LargeLocusCodeAllocation.o \
        $(DIR_OBJ)/Utilities.o \
        $(DIR_OBJ)/DeBruijnGraph.o \
        $(DIR_OBJ)/DeBruijnElement.o \
        $(DIR_OBJ)/basic.o \
        $(DIR_OBJ)/binarykMer.o \
        $(DIR_OBJ)/Hsh.o \
        $(DIR_OBJ)/GraphAligner.o \
        $(DIR_OBJ)/GraphAlignernonAffine.o \
        $(DIR_OBJ)/GraphAlignerAffine.o \
        $(DIR_OBJ)/GraphAlignerendsFree.o \
        $(DIR_OBJ)/GraphAlignerUnique.o \
        $(DIR_OBJ)/GraphAndEdgeIndex.o \
        $(DIR_OBJ)/UniqueAlignerTests.o \
        $(DIR_OBJ)/GraphAndIndex.o \
        $(DIR_OBJ)/AlignerTests.o \
        $(DIR_OBJ)/VirtualNW.o \
        $(DIR_OBJ)/VirtualNWUnique.o
        
#
# list executable file names
#
EXECS = MHC-PRG

#
# compile and link
#
default:
	@echo
	@echo " to build:"
	@echo "    make all"
	@echo
	@echo " to clean:"
	@echo "    make clean"
	@echo "    make realclean"
	@echo

all: $(EXECS)

$(EXECS): $(OBJS)
	$(foreach EX, $(EXECS), $(COMPILE) $(EX).cpp -c -o $(DIR_OBJ)/$(EX).o;)
	$(foreach EX, $(EXECS), $(COMPILE) $(OBJS) $(DIR_OBJ)/$(EX).o -o $(DIR_BIN)/$(EX) $(LIBS);)

$(DIR_OBJ)/%.o: %.cpp %.h
	$(COMPILE) $< -c -o $@


#
# odds and ends
#
clean:
	/bin/rm $(DIR_OBJ)/*

realclean: clean
	/bin/rm $(DIR_BIN)/*
