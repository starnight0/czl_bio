CC=g++
#CXXFLAGS=-DPTHREAD
CXXFLAGS=-I/home/chenz11/include
CXXLIBS=-L/home/chenz11/lib -lz -lpthread 
BAMINC=-I/home/chenz11/include
BAMLIB=-lhts -lbam
CPPFLAGS=-g -fno-inline
CPPLIBS=
# CXXLIBS=-lpthread -lboost_system -lboost_filesystem
# include seq/makefile
obj=czl_common.o czl_bio_base.o czl_io.o czl_bam.o czl_bio_seq.o
sub=seq
libczl_bio: ${obj} $(sub)
	g++ ${obj} seq/bit2_seq.o ${CXXFLAGS} ${CPPFLAGS} $(BAMINC) $(CXXLIBS) -fPIC -shared -o ~/lib/libczl_bio.so
$(obj): %.o: %.cpp
	g++ ${CXXFLAGS} ${CPPFLAGS} -c -fPIC $< -o $@
seq:
	cd seq && $(MAKE)

.PHONY: clean
clean: 
	-rm *.o */*.o
