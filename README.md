# Pairwise Hierarchical Model (Underconstruction)
A Bayesian hierarchical model to map causal interactions between regulatory elements in the genome, that incorporates techniques from Mendelian Randomisation

## Bugs

Non sequential input for hm.c.

## How to build & install
**Please make sure HTSlib and CLAPACK are installed in your environment** (if you don't have them, then see below for installation tips).

To build and install the pairwise hierarchical model, firstly go to the _source_ directory (*src*), then set environment variables appropriately to point to the CLAPACK and HTSLIB.  Finally use "make" to build and install *phm* which will be installed in "$PHMDIR/bin".

        PHMDIR=/path/to/phmdir/
        cd $PHMDIR/src
        # Not run!  Please export your environment.
        export CFLAGS="-I/path/to/your/CLAPACK-*.*.*.*/INCLUDE -I/path/to/your/CLAPACK-*.*.*.*/F2CLIBS"
        export LDFLAGS="-L/path/to/your/CLAPACK-*.*.*.* -L/path/to/your/CLAPACK-*.*.*.*/F2CLIBS"
        make
        make install

## Workflow

![workflow](https://github.com/natsuhiko/Images/blob/master/workflow.png)

## Bayes factor calculation

PHM takes Bayes factors (BFs) of QTL associations as an input data. Here we describe how to compute BFs from normalised read counts and variant information in VCF format.

	bayeslm -g /path/to/your/VCF/file.gz -w 1000000 -j 10 -f /path/to/your/peak.bed.gz

Note that the BED file of the peak annotation has to be tabix indexed. The regional Bayes factors for PHM also calculated with the parameter esitimate from the hierarchical model.

	bayeslm -g /path/to/your/VCF/file.gz -w 1000000 -j 10 -f /path/to/your/peak.bed.gz -p1 variant_level_prior.gz -p2 peak_level_prior.gz

The function provides the regional bayes factors for all peak k (>j) in the cis-window.

## Model fitting

Pairwise hierarchical model uses QTL signal to map causal interaction between regulatory elements (hereafter we refer those elemnts as *peaks*). The model employs 2 stage optimisation: the first step is fitting a standard hierarchical model to estimate the variant-level and peak-level prior probabilities under the assumption that peaks are independent; then the second step is fitting the pairwise hierarchical model using the parameter estimate in the first stage. At each model fitting stage, the Bayes factors of genetic associations are required. The next section provides the information. 

	# First stage - standard hierarchical model
	hm  bf1.gz -v variant_level_prior.gz -f feature_level_prior.gz

	# Second stage - pairwise hierarchical model
	phm bf2.gz -p posterior_prob.gz -c coef.gz

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

You may also need to obtain GSL (GNU Scientific Library) from http://www.gnu.org/software/gsl/ (if it is not installed).  Then, compile it like

        tar zxvf gsl-*.*.tar.gz
        cd gsl-*.*
        ./configure --prefix=$PWD
        make
        make install

