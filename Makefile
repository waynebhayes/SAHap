CC=gcc
CXX=g++
CXXFLAGS = -g -I"include" -Wall -std=c++11 -O3 #-ggdb -pg

sahap: src/main.o src/Allele.o src/Haplotype.o src/Genome.o src/InputReader.o src/utils.o
	g++ -std=c++11 -o sahap src/*.o

src/main.o: src/main.cpp
src/Allele.o: src/Allele.cpp
src/Haplotype.o: src/Haplotype.cpp
src/Genome.o: src/Genome.cpp
src/InputReader.o: src/InputReader.cpp
src/utils.o: src/utils.cpp

clean:
	/bin/rm -f *.o */*.o sahap sahap.exe
