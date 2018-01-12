PROGFLAG = -DNIL -DFULLBF -DPEN
CC=gcc
CFLAGS := $(CFLAGS) -std=gnu99 -fpic -g -O2
LDFLAGS := $(LDFLAGS) -lgsl -llapack -lblas -lm -lz -lhts

PROGRAM = phm
HM = hm
BL = bayeslm
OBJS = $(PROGRAM).o util.o usage.o
OBJS_HM = $(HM).o util.o usage.o
OBJS_BL = $(BL).o util_bayeslm.o usage.o getLogBF.o

.SUFFIXES: .c .o

.PHONY: all
all: $(PROGRAM) $(HM) $(BL)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(PROGFLAG) -o $(PROGRAM) $^

$(HM): $(OBJS_HM)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(HM) $^

$(BL): $(OBJS_BL)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(BL) $^

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
