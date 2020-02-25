CPP       = g++
CPPARGS   = -O3 -m64 -std=c++14 -Wall -Wextra -pedantic
SRCPATH		= ./src/
BINPATH		= ./bin/
DATPATH		= ./dat/
## GRBPATH =  /util/academic/gurobi/gurobi900/linux64
GRBPATH =  /opt/gurobi900/linux64/
INCGRB    = -I$(GRBPATH)/include/
CPPLIBGRB = -L$(GRBPATH)/lib/ -lgurobi_c++ -lgurobi90 $(CPPSTDLIB) -lpthread -lm


all:
	$(CPP) $(CPPARGS) $(SRCPATH)DataHandler.cpp \
	$(SRCPATH)PartitionScheme.cpp \
	$(SRCPATH)STEFormulation.cpp \
	$(SRCPATH)BendersCuts.cpp \
	$(SRCPATH)DualFormulation.cpp \
	$(SRCPATH)GlobalMC.cpp \
	$(SRCPATH)SuperCutFormulation.cpp \
	$(SRCPATH)SubtourCuts.cpp \
	-o $(BINPATH)main $(SRCPATH)main.cpp $(INCGRB) $(CPPLIBGRB)

clean:
	rm -rf $(BINPATH)*.o $(BINPATH)*.dSYM $(BINPATH)main

run:
	$(BINPATH)main
