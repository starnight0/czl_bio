CC=g++
#CXXFLAGS=-DPTHREAD
CXXFLAGS=-I/home/chenz11/include
CXXLIBS=-L/home/chenz11/lib -lz -lpthread 
BAMINC=-I/home/chenz11/include
BAMLIB=-lhts -lbam
CPPFLAGS=-g
CPPLIBS=
# CXXLIBS=-lpthread -lboost_system -lboost_filesystem
obj=czl_common.o czl_bio_base.o czl_io.o czl_bam.o
libczl_bio: ${obj}
	g++ ${obj} ${CXXFLAGS} ${CPPFLAGS} $(BAMINC) $(CXXLIBS) -fPIC -shared -o ~/lib/libczl_bio.so
czl_common.o: czl_common.cpp
	g++ ${CXXFLAGS} ${CPPFLAGS} -c -fPIC czl_common.cpp
czl_bio_base.o: czl_bio_base.cpp
	g++ ${CXXFLAGS} ${CPPFLAGS} -c -fPIC czl_bio_base.cpp
czl_io.o: czl_io.cpp
	g++ ${CXXFLAGS} ${CPPFLAGS} -c -fPIC czl_io.cpp
czl_bam.o: czl_bam.cpp
	g++ ${CXXFLAGS} ${CPPFLAGS} -c -fPIC czl_bam.cpp

.PHONY: clean
clean: 
	rm *.o
