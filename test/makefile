CXXFLAGS=-L$(HOME)/include 
CXXLIBS=-L$(HOME)/lib -lczl_bio -pthread -lbam -lhts
CPPFLAGS=-g -fno-inline-functions
CPPLIBS=
test: test.cc
	g++ -g -fno-inline-functions -o $@ $^ -lczl_bio2 -lboost_system -lgtest -pthread
test.StringUtility: test.StringUtility.cpp
	g++ $^ -g -fno-inline-functions $(LIBS) $(CXXFLAGS) $(CPPFLAGS) -o $@ -lczl_bio -lgtest -pthread
test.seq: test.seq.cpp
	g++ $^ $(CXXFLAGS) $(CPPFLAGS) -o $@ $(CXXLIBS) $(CPPLIBS)
clear:
	rm test
