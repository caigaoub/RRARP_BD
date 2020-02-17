CPP       = g++
CPPARGS   = -O3 -m64 -std=c++1y -Wall -Wextra -pedantic
SRCPATH		= ./src/
BINPATH		= ./bin/
DATPATH		= ./dat/
# GRBPATH   = /opt/gurobi900/linux64
GRBPATH =  /util/academic/gurobi/gurobi900/linux64
INCGRB    = -I$(GRBPATH)/include/
CPPLIBGRB = -L$(GRBPATH)/lib/ -lgurobi_c++ -lgurobi90 $(CPPSTDLIB) -lpthread -lm
# INSTANCE	= 7 $(DATPATH)sample_n_25.txt

all:
	$(CPP) $(CPPARGS) $(SRCPATH)DataHandler.cpp \
	$(SRCPATH)PartitionScheme.cpp \
	$(SRCPATH)STEFormulation.cpp \
	$(SRCPATH)BendersCuts.cpp \
	$(SRCPATH)DualFormulation.cpp \
	$(SRCPATH)GlobalMC.cpp \
	-o $(BINPATH)main $(SRCPATH)main.cpp $(INCGRB) $(CPPLIBGRB)

clean:
	rm -rf $(BINPATH)*.o $(BINPATH)*.dSYM $(BINPATH)main

run:
	$(BINPATH)main
