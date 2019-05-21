CPP       = g++
CPPARGS   = -O3 -m64 -std=c++14 -Wall -Wextra -pedantic
SRCPATH		= ./
BINPATH		= ./
DATPATH		= ./
GRBPATH   = /opt/gurobi801/linux64
INCGRB    = -I$(GRBPATH)/include/
CPPLIBGRB = -L$(GRBPATH)/lib/ -lgurobi_c++ -lgurobi80 $(CPPSTDLIB) -lpthread -lm

all:
	$(CPP) $(CPPARGS) $(SRCPATH)InstanceGenerator.cpp -g -o $(BINPATH)main $(SRCPATH)main.cpp $(INCGRB) $(CPPLIBGRB)

clean:
	rm -rf $(BINPATH)*.o $(BINPATH)*.dSYM $(BINPATH)main
