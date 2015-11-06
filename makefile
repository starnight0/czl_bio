CC=g++
CXXFLAGS=-DPTHREAD
CXXLIBS=-lpthead
obj=czl_common.o czl_bio_base.o czl_io.o
libczl_bio2: ${obj}
	g++ -g ${obj} ${CXXFLAGS} -fPIC -shared -o /usr/local/lib/libczl_bio2.so
czl_common.o: czl_common.h czl_common.cc
	g++ -g ${CXXFLAGS} -c -fPIC czl_common.cc
czl_io.o: czl_io.h czl_io.cc czl_common.h czl_bio_base.h
	g++ -g ${CXXFLAGS} -I/home/chzelin/bioprogram/samtools-0.1.19/include -c -fPIC czl_io.cc
czl_bio_annot.o: czl_bio_annot.h czl_bio_annot.cc
	g++ -g ${CXXFLAGS} -c -fPIC czl_bio_annot.cc
czl_log.o: czl_log.h czl_log.cc
	g++ -g ${CXXFLAGS} -c -fPIC czl_log.cc -lpthread
czl_blast.o: czl_blast.h czl_blast.cc
	g++ -g ${CXXFLAGS} -c -fPIC czl_blast.cc -lpthread
czl_bio_base.o: czl_bio_base.h czl_bio_base.cc czl_common.o
	g++ -g ${CXXFLAGS} -c -fPIC czl_bio_base.cc
czl_bam.o: czl_bam.h czl_bam.cc
	g++ -g ${CXXFLAGS} -c -fPIC czl_bam.cc czl_bam.h ${CXXLIBS}
