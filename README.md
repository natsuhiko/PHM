# Pairwise Hierarchical Model
A Bayesian hierarchical model to map causal interactions between regulatory elements in the genome, that incorporates techniques from Mendelian Randomisation

## How to build & install
**Please make sure HTSLIB and CLAPACK are installed in your environment** (if you don't have them, then see below for installation tips).

To build and install the pairwise hierarchical model, firstly go to the _source_ directory (*src*), then set environment variables appropriately to point to the CLAPACK and HTSLIB.  Finally use "make" to build and install *phm* which will be installed in "$PHMDIR/bin".

        PHMDIR=/path/to/phmdir/
        cd $PHMDIR/src
        # Not run!  Please export your environment.
        export CFLAGS="-I/path/to/your/CLAPACK-*.*.*.*/INCLUDE -I/path/to/your/CLAPACK-*.*.*.*/F2CLIBS"
        export LDFLAGS="-L/path/to/your/CLAPACK-*.*.*.* -L/path/to/your/CLAPACK-*.*.*.*/F2CLIBS"
        make
        make install

## Model fitting

Pairwise hierarchical model uses QTL signal to map causal interaction between regulatory elements (hereafter we refer those elemnts as *peaks*). The model employs 2 stage optimisation: the first step is fitting a standard hierarchical model to estimate the variant-level and peak-level prior probabilities under the assumption that peaks are independent; then the second step is fitting the pairwise hierarchical model using the parameter estimate in the first stage. At each model fitting stage, the Bayes factors of genetic associations are required. The next section provides the information. 

	hm  bf1.gz -v variant_level.gz -f peak_level.gz
	phm bf2.gz -p posterior_prob.gz -c coef.gz

## Bayes factor calculation

PHM takes Bayes factors of QTL associations.

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
