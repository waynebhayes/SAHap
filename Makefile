CC=gcc
CXX=g++
CXXFLAGS = -I"include" -Wall -std=c++11 -O3 #-ggdb -pg

sahap: main.o src/Allele.o src/Haplotype.o src/Genome.o src/InputReader.o src/utils.o
	g++ -std=c++11 -o sahap *.o src/*.o

main.o: main.cpp
src/Allele.o: src/Allele.cpp
src/Haplotype.o: src/Haplotype.cpp
src/Genome.o: src/Genome.cpp
src/InputReader.o: src/InputReader.cpp
src/utils.o: src/utils.cpp

clean:
	/bin/rm -f *.o */*.o
