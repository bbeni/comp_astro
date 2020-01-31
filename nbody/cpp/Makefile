CC=g++
CFLAGS=-I. -std=c++14 -O3 -ffast-math

all: hermite.o
	$(CC) $(CFLAGS) -o hermite hermite.o

hermite.o: hermite.cpp
	$(CC) $(CFLAGS) -c hermite.cpp

clean:
	rm hermite hermite.o 

rmsnaps:
	rm snapshot*.csv