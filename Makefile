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
CXXFLAGS := $(CXXFLAGS) '-DOBJECTIVE=OBJ_$(OBJECTIVE)'
INCLUDES := include/Allele.hpp include/Genome.hpp include/Haplotype.hpp include/InputReader.hpp include/ScoringModel.hpp include/types.hpp include/utils.hpp

sahap.$(OBJECTIVE): src/main.o src/Allele.o src/Haplotype.o src/Genome.o src/InputReader.o src/utils.o
	g++ -std=c++11 -o sahap.$(OBJECTIVE) src/*.o

MEC:
	grep -q MEC .last-made || rm -f *.o */*.o && echo MEC > .last-made
	$(MAKE) 'OBJECTIVE=MEC'

Poisson:
	grep -q Poisson .last-made || rm -f *.o */*.o && echo Poisson > .last-made
	$(MAKE) 'OBJECTIVE=Poisson'

src/main.o: src/main.cpp $(INCLUDES)
src/Allele.o: src/Allele.cpp $(INCLUDES)
src/Haplotype.o: src/Haplotype.cpp $(INCLUDES)
src/Genome.o: src/Genome.cpp $(INCLUDES)
src/InputReader.o: src/InputReader.cpp $(INCLUDES)
src/utils.o: src/utils.cpp $(INCLUDES)

parallel: src/parallel.c
	gcc -o parallel src/parallel.c

clean:
	/bin/rm -f parallel sahap.MEC sahap.Poisson *.o */*.o *.exe
