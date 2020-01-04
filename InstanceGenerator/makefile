CPP       = g++
CPPARGS   = -O3 -m64 -std=c++14 -Wall -Wextra -pedantic
SRCPATH		= ./src/
BINPATH		= ./bin/
DATPATH		= ./
GRBPATH   = /opt/gurobi900/linux64
INCGRB    = -I$(GRBPATH)/include/
INCBOOST  = -I /usr/include/
CPPLIBGRB = -L$(GRBPATH)/lib/ -lgurobi_c++ -lgurobi90  $(CPPSTDLIB) -lpthread -lm
CPPLIBBOOST = -L /usr/lib/  -lboost_system -lboost_filesystem


all:
	$(CPP) $(CPPARGS) $(SRCPATH)InstanceGenerator.cpp $(SRCPATH)Reader.cpp -g -o $(BINPATH)main  $(SRCPATH)main.cpp $(INCGRB) $(INCBOOST) $(CPPLIBGRB) $(CPPLIBBOOST)

# clean:
#	rm -rf $(BINPATH)*.o $(BINPATH)*.dSYM $(BINPATH)main
