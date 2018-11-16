PROGFLAG = -DNIL -DFULLBF -DPEN
CC=gcc
CFLAGS := $(CFLAGS) -std=gnu99 -I/usr/include -I/usr/include/gsl -fpic -g -O2
LDFLAGS := $(LDFLAGS) -lgslcblas -lgsl -llapack -lblas  -lm -lz -lpthread -lf2c -lhts

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
	$(CC) $(CFLAGS) $(PROGFLAG) -o $(PROGRAM) $^ $(LDFLAGS)

$(HM): $(OBJS_HM)
	$(CC) $(CFLAGS) -o $(HM) $^ $(LDFLAGS)

$(BL): $(OBJS_BL)
	$(CC) $(CFLAGS) -o $(BL) $^ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) -c $< $(PROGFLAG)

usage.c:
	sh makeUsage.sh


install: $(PROGRAM) $(HM) $(BL)
	mkdir -p ../bin
	mv phm ../bin
	mv hm ../bin
	mv bayeslm ../bin

clean:
	rm *o
	rm usage.c
