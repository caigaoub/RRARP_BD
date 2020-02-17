CPP             = g++
CPPARGS         = -o3 -m64 -g -std=c++11
SRCPATH		= ./src/
BINPATH		= ./bin/
GRBPATH         = /util/academic/gurobi/gurobi650/linux64
INCGRB          = $(GRBPATH)/include/
INCBOOST	= /projects/academic/josewalt/ningjwe/cpplib/boost_1_60_0/
CPPLIBGRB       = -L$(GRBPATH)/lib/ -lgurobi_c++ -lgurobi65

all: 
	$(CPP) $(CPPARGS) $(SRCPATH)Callback_supercuts.cpp $(SRCPATH)DataReader.cpp $(SRCPATH)Graph.cpp $(SRCPATH)Graph_MinC.cpp $(SRCPATH)Graph_Sub.cpp \
	$(SRCPATH)MainTree.cpp $(SRCPATH)MF_Model.cpp $(SRCPATH)Model_Coll.cpp $(SRCPATH)Model_Sol.cpp $(SRCPATH)model_supercuts.cpp \
	$(SRCPATH)MSTBlockerModel.cpp $(SRCPATH)MSTCuts.cpp $(SRCPATH)Node.cpp $(SRCPATH)ST.cpp $(SRCPATH)StaticParam.cpp $(SRCPATH)StrManip.cpp $(SRCPATH)SysInfo.cpp \
	$(SRCPATH)TCLimits.cpp $(SRCPATH)Timer.cpp -o $(BINPATH)main $(SRCPATH)main.cpp -I$(INCGRB) -I$(INCBOOST) $(CPPLIBGRB) -lpthread -lm 

clean:
	rm -rf $(BINPATH)*.o $(BINPATH)*.dSYM $(BINPATH)main
        
