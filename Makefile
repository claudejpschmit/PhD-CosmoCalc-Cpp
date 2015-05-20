#################################################
# This file should just make a bunch of objects,#
# to simulate what class does.                  #
#################################################

MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build
LIBRARIES = lib/
.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base

vpath %.cpp src:main:$(LIBRARIES)ALGLIB_src
vpath %.c $(LIBRARIES)CLASS_source:$(LIBRARIES)CLASS_tools:$(LIBRARIES)CLASS_main
vpath %.o build
vpath .base build

# Compiler required for c++ code.
CXX = g++ -Wall -std=c++11
# Compiler required for c code.
CC = gcc -Wall

OPTFLAG = -O4

# This decreases the precision of the bessel functions significantly if 
# added to the compilation of the files containing boost->sph_bess(l,x).
OPTFLAG_CLASS = -ffast-math
OMPFLAG = -fopenmp
CCFLAG = -g -fPIC
LDFLAG = -g -fPIC

# leave blank to compile without HyRec, or put path to HyRec directory
# (with no slash at the end: e.g. CLASS_hyrec or ../CLASS_hyrec)
# TODO: This is not working yet!!!
HYREC = $(LIBRARIES)CLASS_hyrec 
# automatically add external programs if needed. First, initialize to blank.
EXTERNAL =

CCFLAG += -D__CLASSDIR__='"$(MDIR)"'
INCLUDES = -I../include
INCLUDES += -I../$(LIBRARIES)CLASS_include
INCLUDES += -I../$(LIBRARIES)ALGLIB_include

# These lines seem to be unnecessary, but I leave them in anyways.
INCLUDES += -I/usr/include/boost
LINKER = -L/usr/include/boost #-lboost_filesystem

# eventually update flags for including HyRec
ifneq ($(HYREC),)
vpath %.c $(HYREC)
CCFLAG += -DHYREC
#LDFLAGS += -DHYREC
INCLUDES += -I../$(LIBRARIES)CLASS_hyrec
EXTERNAL += hyrectools.o helium.o hydrogen.o history.o
endif

%.o: %.cpp .base
	cd $(WRKDIR);$(CXX) $(OPTFLAG) $(OMPFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

# This line creates the CLASS objects.
%.o: %.c .base
	cd $(WRKDIR);$(CC) $(OPTFLAG) $(OPTFLAG_CLASS) $(OMPFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

ALGLIB = alglibinternal.o alglibmisc.o ap.o dataanalysis.o diffequations.o fasttransforms.o integration.o interpolation.o linalg.o optimization.o solvers.o specialfunctions.o statistics.o
TOOLS = growTable.o dei_rkck.o sparse.o evolver_rkck.o  evolver_ndf15.o arrays.o parser.o quadrature.o hyperspherical.o common.o
SOURCE = input.o background.o thermodynamics.o perturbations.o primordial.o nonlinear.o transfer.o spectra.o lensing.o
CLASS = class.o
OUTPUT = output.o

SRC = CosmoBasis.o CosmologyCalculatorClass.o Engine.o ClassEngine.o
MAIN = Main.o

all: calc class_test

calc: $(SRC) $(SOURCE) $(TOOLS) $(OUTPUT) $(ALGLIB) $(EXTERNAL) $(MAIN) 
	cd $(MDIR);$(CXX) $(OPTFLAG) $(OPTFLAG_CLASS) $(OMPFLAG) $(LDFLAG) $(LINKER) -o calc $(addprefix build/, $(notdir $^)) -lm

class_test: $(SOURCE) $(TOOLS) $(OUTPUT) $(EXTERNAL) $(CLASS) 
	cd $(MDIR);$(CC) $(OPTFLAG) $(OPTFLAG_CLASS) $(OMPFLAG) $(LDFLAG) $(LINKER) -o class_test $(addprefix build/, $(notdir $^)) -lm


clean: .base
	rm -rf $(WRKDIR);
	rm -f libclass.a
