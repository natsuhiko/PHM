CLAPACKINC=/nfs/users/nfs_n/nk5/Applications/CLAPACK-3.1.1.1/INCLUDE
CLAPACKLIBS=/nfs/users/nfs_n/nk5/Applications/CLAPACK-3.1.1.1
F2CLIBS=/nfs/users/nfs_n/nk5/Applications/CLAPACK-3.1.1.1/F2CLIBS
GSLLIBS=/nfs/users/nfs_n/nk5/Applications/gsl-1.16/lib
GSLINC=/nfs/users/nfs_n/nk5/Applications/gsl-1.16/include/gsl
HTSLIB=/nfs/users/nfs_n/nk5/Applications/htslib-1.3

PROGFLAG = -DFIT # -DPROOT
CC=gcc
CFLAGS := -std=gnu99 -I$(CLAPACKINC) -I$(GSLINC) -I$(HTSLIB) -fpic -g -O2
LDFLAGS := -L$(CLAPACKLIBS) -L$(GSLLIBS) -L$(HTSLIB) -lhts -llapack -lblas -lm -lz -lgsl

PROGRAM = hierarchical2

OBJS = $(PROGRAM).o util.o

.SUFFIXES: .c .o

.PHONY: all
all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) -o $(PROGRAM) $^ $(LDFLAGS) $(PROGFLAG)

.c.o:
	$(CC) $(CFLAGS) -c $< $(PROGFLAG)

clean:
	rm *o
