LIBS=-L/home/chzelin/lib
CPPFLAGS=-I~/include
test: test.cc
	g++ -g -fno-inline-functions -o $@ $^ -lczl_bio2 -lboost_system -lgtest -pthread
test.StringUtility: test.StringUtility.cpp
	g++ $^ -g -fno-inline-functions $(LIBS) $(CXXFLAGS) $(CPPFLAGS) -o $@ -lczl_bio -lgtest -pthread
clear:
	rm test
