CC=gcc
CXX=g++
CXXFLAGS = -g -I"include" -Wall -std=c++11 -O3 #-ggdb -pg

ifeq ($(OBJECTIVE),MEC)
    found=1
else ifeq ($(OBJECTIVE),Poisson)
    found=1
else
    OBJECTIVE=MEC
endif
CXXFLAGS := $(CXXFLAGS) '-DOBJECTIVE=$(OBJECTIVE)'

sahap.$(OBJECTIVE): src/main.o src/Allele.o src/Haplotype.o src/Genome.o src/InputReader.o src/utils.o
	g++ -std=c++11 -o sahap.$(OBJECTIVE) src/*.o

MEC:
	make clean
	$(MAKE) 'OBJECTIVE=MEC'

Poisson:
	make clean
	$(MAKE) 'OBJECTIVE=Poisson'

src/main.o: src/main.cpp
src/Allele.o: src/Allele.cpp
src/Haplotype.o: src/Haplotype.cpp
src/Genome.o: src/Genome.cpp
src/InputReader.o: src/InputReader.cpp
src/utils.o: src/utils.cpp

clean:
	/bin/rm -f *.o */*.o
