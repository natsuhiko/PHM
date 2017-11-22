CLAPACKINC=/nfs/users/nfs_n/nk5/Applications/CLAPACK_GCC/CLAPACK-3.1.1.1/INCLUDE
CLAPACKLIBS=/nfs/users/nfs_n/nk5/Applications/CLAPACK_GCC/CLAPACK-3.1.1.1
F2CLIBS=/nfs/users/nfs_n/nk5/Applications/CLAPACK_GCC/CLAPACK-3.1.1.1/F2CLIBS
GSLLIBS=/nfs/users/nfs_n/nk5/Applications/gsl-1.16/lib
GSLINC=/nfs/users/nfs_n/nk5/Applications/gsl-1.16/include/gsl
HTSLIB=/nfs/users/nfs_n/nk5/Applications/htslib-1.3

PROGFLAG = -DFIT # -DPROOT
CC=gcc
CFLAGS := $(CFLAGS) -std=gnu99 -fpic -g -O2
LDFLAGS := $(LDFLAGS) -llapack -lblas -lm -lz -lhts

PROGRAM = hierarchical2

OBJS = $(PROGRAM).o util.o

.SUFFIXES: .c .o

.PHONY: all
all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(PROGFLAG) -o $(PROGRAM) $^

.c.o:
	$(CC) $(CFLAGS) -c $< $(PROGFLAG)

clean:
	rm *o
