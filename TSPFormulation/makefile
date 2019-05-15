CPP       = g++
CPPARGS   = -O3 -m64 -std=c++14 -Wall -Wextra -pedantic
SRCPATH		= ./src/
BINPATH		= ./bin/
DATPATH		= ./dat/
GRBPATH   = /opt/gurobi801/linux64
INCGRB    = -I$(GRBPATH)/include/
CPPLIBGRB = -L$(GRBPATH)/lib/ -lgurobi_c++ -lgurobi80 $(CPPSTDLIB) -lpthread -lm

all:
	$(CPP) $(CPPARGS) $(SRCPATH)SubtourElimCuts.cpp \
	$(SRCPATH)STEFormulation.cpp \
	$(SRCPATH)costInfo.cpp -o $(BINPATH)MainRun $(SRCPATH)main.cpp $(INCGRB) $(CPPLIBGRB)

clean:
	rm -rf $(BINPATH)*.o $(BINPATH)*.dSYM $(BINPATH)main
