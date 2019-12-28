CC=gcc
CXX=g++
CXXFLAGS = -I"include" -Wall -std=c++11 -O3 #-ggdb -pg

sahap: main.o src/Allele.o src/Chromosome.o src/Genome.o src/InputReader.o
	g++ -std=c++11 -o sahap *.o src/*.o

main.o: main.cpp
src/Allele.o: src/Allele.cpp
src/Chromosome.o: src/Chromosome.cpp
src/Genome.o: src/Genome.cpp
src/InputReader.o: src/InputReader.cpp

clean:
	/bin/rm -f *.o */*.o
