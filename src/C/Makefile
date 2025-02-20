CC=gcc -O2 #-ggdb

sahap-c: bitvec.o sahap-c.o misc.o sets.o
	$(CC) -o sahap-c *.o -lm

clean:
	rm -f *.o *.exe*
