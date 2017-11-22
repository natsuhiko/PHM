# Pairwise Hierarchical Model
A Bayesian hierarchical model that incorporates techniques from Mendelian Randomisation

## How to build & install
**Please make sure HTSLIB, CLAPACK and GSL (GNU Scientific Library) are installed in your environment** (if you don't have them, then see below for installation tips).

## Installation tips for CLAPACK and GSL

You first need to get the latest library from http://www.netlib.org/clapack/.  Then, compile it like

        tar zxvf clapack.tgz
        cd CLAPACK-3.2.1
        mv make.inc.example make.inc
        make

When it has been done, you will find three archive files in the directory which need to be either linked or renamed, such as

        ln -s lapack_LINUX.a liblapack.a
        ln -s tmglib_LINUX.a libtmglib.a
        ln -s blas_LINUX.a libblas.a

before compiling the pairwise hierarchical model.

You may also need to get GSL (GNU Scientific Library) from http://www.gnu.org/software/gsl/ (if it is not installed).  Then, compile it like

        tar zxvf gsl-*.*.tar.gz
        cd gsl-*.*
        ./configure --prefix=$PWD
        make
        make install
