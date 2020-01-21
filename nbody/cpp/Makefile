CC=g++ -O3 -ffast-math
CFLAGS=-I.

all: first_task.o
	$(CC) -o first_task first_task.o

first_task.o: first_task.cpp
	$(CC) -c first_task.cpp

clean:
	rm first_task first_task.o 