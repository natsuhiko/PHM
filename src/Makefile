PROGFLAG = -DNIL -DFULLBF
CC=gcc
CFLAGS := $(CFLAGS) -std=gnu99 -fpic -g -O2
LDFLAGS := $(LDFLAGS) -llapack -lblas -lm -lz -lhts

PROGRAM = phm
OBJS = $(PROGRAM).o util.o usage.o

.SUFFIXES: .c .o

.PHONY: all
all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(PROGFLAG) -o $(PROGRAM) $^

.c.o:
	$(CC) $(CFLAGS) -c $< $(PROGFLAG)

usage.c:
	sh makeUsage.sh


install:
	mkdir -p ../bin
	mv phm ../bin

clean:
	rm *o
	rm usage.c
