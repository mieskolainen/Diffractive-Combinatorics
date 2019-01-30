# Makefile for Diffractive Cross Section analysis
#
#
# Compiling:
#
#        make dictionary (REMEMBER to do this if the ROOT class has been modified)
#		 make
#
# Cleaning:
#
#        make clean
#        make superclean
#
# mikael.mieskolainen@cern.ch, 2018
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
#
# USE [TABS] for intendation while modifying this!


# -----------------------------------------------------------------------
# Basic paths

# HepMC2 installation path
HEPMC2DIR            = /home/user/cernbox/Granite/HepMC-2.06.09/install

# ROOT installation path
ROOTLIBDIR           = $(shell root-config --libdir)

# AliROOT/AliPhysics path
ALIROOT              = /home/user/alice/sw/ubuntu1604_x86-64/AliRoot/0-1
ALIPHYSICS           = /home/user/alice/sw/ubuntu1604_x86-64/AliPhysics/0-1

# RooUnfold local path
ROOUNFOLD            = /home/user/alice/RooUnfold

# ------------------------------------------------------------------------
# Libraries (compiled)

# HepMC2 libraries
HEPMC2LIB            = -L$(HEPMC2DIR)/lib -lHepMC
HEPMC2FIOLIB         = -L$(HEPMC2DIR)/lib -lHepMCfio

# ROOT libraries
ROOTlib              = -L$(ROOTLIBDIR) -lCore -lCint -lRIO -lNet \
                       -lHist -lGraf -lGraf3d -lGpad -lTree -lRint \
                       -lPostscript -lMatrix -lPhysics -lMathCore \
                       -lThread -lGui -lRooFit -lMinuit -lEG

# ALICE libraries
ALIROOTlib           = -L$(ALIROOT)/lib -lESD -lSTEER -lSTEERBase -lANALYSISalice -lANALYSIS
ALIPHYSICSlib        = -L$(ALIPHYSICS)/lib -lOADB

ROOUNFOLDlib         = -L/home/user/alice/sw/BUILD/RooUnfold-latest/RooUnfold -lRooUnfold


# C++ standard libraries
STANDARDlib          = -pthread -rdynamic -lm -ldl

# Boost libraries
BOOSTLIBS    	     = -lboost_iostreams -lboost_system


# ------------------------------------------------------------------------
# Source and objects and dependencies

SRC_DIR = src
OBJ_DIR = obj

SRC     = $(wildcard $(SRC_DIR)/*.cc)
OBJ     = $(SRC:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)
DEP     = $(OBJ:$(OBJ_DIR)/%.o=.d)


# ------------------------------------------------------------------------
# Header files

INCLUDES   =  -I.
INCLUDES   += -Iinclude
INCLUDES   += -Ilibs
INCLUDES   +=  $(STANDARDlib)
INCLUDES   += -I$(ROOTSYS)/include -I$(ALIROOT)/include -I$(ALIPHYSICS)/include -I$(ROOUNFOLD)/include
INCLUDES   += -I$(HEPMC2DIR)/include


# -----------------------------------------------------------------------------
# External libraries to be linked

LINK_LIBS  =  $(ROOTlib)
LINK_LIBS  += $(ALIROOTlib)
LINK_LIBS  += $(ALIPHYSICSlib)
LINK_LIBS  += $(ROOUNFOLDlib)
LINK_LIBS  += $(HEPMC2LIB)


# ------------------------------------------------------------------------
# Compiler options

CXX       = g++
CXXFLAGS  = -std=c++14 -Wall -fPIC -pipe -O3 -march=native -ftree-vectorize $(INCLUDES)

# Automatic dependencies on
CXXFLAGS += -MMD -MP

# -Wall,           compiler warnings full on
# -march=native,   CPU spesific instruction set usage
# -free-vectorize, Autovectorization on
# -ffast-math,     Heavy floating point optimization
#  (fast but breaks FLOP IEEE standards, do not use!)
# -fPIC,           position independent code (PIC) for shared libraries
# -O2, -O3,        for optimization
# -pipe,           Faster compilation using pipes
# -pg,             Profiling and debugging (MAKES SIGNIFICANT PERFORMANCE HIT)
#
#
# Check your CPU instruction set with: cat /proc/cpuinfo


# -----------------------------------------------------------------------------
.SUFFIXES:      .o .cc
all:	libraries analysis

# Our object files
libraries: $(OBJ)

# Programs
analysis: analysis.o $(OBJ)
	$(CXX) $(CXXFLAGS) $@.o $(OBJ) $(LINK_LIBS) -o $@

# ROOT Dictionary creation
ALIFILE = AliAnalysisTaskDiffCrossSectionsMM.h
dictionary:
	cp ./include/Linkdef.h Linkdef.h
	cp ./include/$(ALIFILE) .
	rootcint -f Dictionary.cc -c -I$(ALIROOT)/include -I$(ALIPHYSICS)/include -L$(ALIROOT)/lib -L$(ALIPHYSICS)/lib $(ALIFILE) Linkdef.h
	cp Dictionary.h ./include/Dictionary.h
	cp Dictionary.cc ./src/Dictionary.cc
	rm Dictionary.* Linkdef.h $(ALIFILE)


# -----------------------------------------------------------------------------
# Compile .o (object) files from a .cc (source) files

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	@echo " "
	@echo "Generating dependencies and compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@


# -----------------------------------------------------------------------------
# Clean up of object and dependency files
clean:
	rm *.o
	rm $(OBJ_DIR)/*.o
	rm *.d
	rm $(OBJ_DIR)/*.d


# -----------------------------------------------------------------------------
# Clean + clean executables + clear ROOT dictionaries
superclean:
	$(MAKE) clean --no-print-directory
	rm ./include/Dictionary.h
	rm ./src/Dictionary.cc


# -----------------------------------------------------------------------------
# Dependencies listed here
-include $(DEP)
