CPP       	= g++
CPPARGS   	= -O3 -m64 -std=c++1y -Wall -Wextra -pedantic
SRCPATH		= ./src/
BINPATH		= ./bin/
DATPATH		= ./dat/
# GRBPATH 	=  /util/academic/gurobi/gurobi900/linux64/
GRBPATH 	=  /opt/gurobi902/linux64/
INCGRB    	= -I$(GRBPATH)/include/
CPPLIBGRB 	= -L$(GRBPATH)/lib/ -lgurobi_c++ -lgurobi90 $(CPPSTDLIB) -lpthread -lm
INCBOOST  	= -I /usr/include/
CPPLIBBOOST 	= -L /usr/lib/  -lboost_system -lboost_filesystem

all:
	$(CPP) $(CPPARGS) $(SRCPATH)DataHandler.cpp \
	$(SRCPATH)PartitionScheme.cpp \
	$(SRCPATH)STEFormulation.cpp \
	$(SRCPATH)BendersCuts.cpp \
	$(SRCPATH)DualFormulation.cpp \
	$(SRCPATH)GlobalMC.cpp \
	$(SRCPATH)SuperCutFormulation.cpp \
	$(SRCPATH)SubtourCuts.cpp \
	$(SRCPATH)TSPModels.cpp \
	-o $(BINPATH)main $(SRCPATH)main.cpp $(INCGRB) $(INCBOOST) $(CPPLIBGRB) $(CPPLIBBOOST)

clean:
	rm -rf $(BINPATH)*.o $(BINPATH)*.dSYM $(BINPATH)main

run:
	$(BINPATH)main
