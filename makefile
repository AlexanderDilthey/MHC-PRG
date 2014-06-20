# LIBRARY SETTINGS - SET AS NECESSARY
# 
# The library settings typically require some tinkering - for reasons beyond me, sometimes one has to include
# the shared object files (.so), and sometimes the .a files (particularly for bamtools).
# Also, curiously, sometimes bamtools requires the explicit inclusion of libz (either as 
# file or just via -lz)
# The following values work for me (see below for an alternative):
#
LIB_BOOST = /home/dilthey/PnP/libs/boost_1_52_0/
INCS = -I$(LIB_BOOST) -IGraph -I/home/dilthey/bamtools/bamtools/include -I/home/dilthey/bamtools/bamtools/src
LIBS = /home/dilthey/PnP/libs/boost_1_52_0/stage/lib/libboost_random.so /home/dilthey/PnP/libs/boost_1_52_0/stage/lib/libboost_filesystem.so /home/dilthey/PnP/libs/boost_1_52_0/stage/lib/libboost_system.so  /home/dilthey/bamtools/bamtools/lib/libbamtools.a /home/dilthey/bamtools/bamtools/lib/libbamtools-utils.so  /home/dilthey/bamtools/zlib-1.2.7/libz.a

# an alternative line (courtesy Peter Humburg, not working for me but for him) is
# LIBS = /home/dilthey/PnP/libs/boost_1_52_0/stage/lib/libboost_random.so /home/dilthey/PnP/libs/boost_1_52_0/stage/lib/libboost_filesystem.so /home/dilthey/PnP/libs/boost_1_52_0/stage/lib/libboost_system.so /home/dilthey/bamtools/bamtools/lib/libbamtools.so /home/dilthey/bamtools/bamtools/lib/libbamtools-utils.a -lz

MKDIR_P = mkdir -p

.PHONY: directories
	
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

VPATH = Data:Graph:NextGen:hash:hash/deBruijn:hash/sequence:GraphAligner:GraphAlignerUnique:readFilter
        
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
        $(DIR_OBJ)/readSimulator.o \
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
        $(DIR_OBJ)/readFilter.o \
        $(DIR_OBJ)/GraphAndIndex.o \
        $(DIR_OBJ)/AlignerTests.o \
        $(DIR_OBJ)/VirtualNW.o \
        $(DIR_OBJ)/VirtualNWUnique.o
        
#
# list executable file names
#
EXECS = MHC-PRG

OUT_DIR = ../obj ../bin

directories: ${OUT_DIR}


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

all: directories $(EXECS)

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

${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}

