CPP       = g++
CPPARGS   = -O3 -m64 -std=c++1y -Wall -Wextra -pedantic
SRCPATH		= ./src/
BINPATH		= ./bin/
DATPATH		= ./dat/
GRBPATH =  /util/academic/gurobi/gurobi900/linux64
INCGRB    = -I$(GRBPATH)/include/
CPPLIBGRB = -L$(GRBPATH)/lib/ -lgurobi90 -lgurobi_c++ $(CPPSTDLIB) -lpthread -lm


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
