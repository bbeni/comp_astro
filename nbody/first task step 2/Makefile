CC=g++
CFLAGS=-I. -std=c++14 -O3 -ffast-math

all: first_task.o
	$(CC) $(CFLAGS) -o first_task first_task.o

first_task.o: first_task.cpp
	$(CC) $(CFLAGS) -c first_task.cpp

clean:
	rm first_task first_task.o 