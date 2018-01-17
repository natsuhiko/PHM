# Pairwise Hierarchical Model (Underconstruction)
A Bayesian hierarchical model to map causal interactions between regulatory elements in the genome, that incorporates techniques from Mendelian Randomisation.

## How to build & install
**Please make sure GSL, HTSlib and CLAPACK are installed in your environment** (if you don't have them, then see below for installation tips).

To build and install the pairwise hierarchical model, firstly go to the _source_ directory (*src*), then set environment variables appropriately to point to the GSL, CLAPACK and HTSLIB.  Finally use "make" to build and install *phm* which will be installed in "$PHMDIR/bin".

	PHMDIR=/path/to/phmdir/
	cd $PHMDIR/src
	# Not run!  Please export your environment.
	export CFLAGS="-I/usr/include -I/path/to/your/htslib-1.* -I/path/to/your/CLAPACK-3.*.*.*/INCLUDE"
	export LDFLAGS="-L/usr/lib -L/path/to/your/htslib-1.* -L/path/to/your/CLAPACK-3.*.*.*"
	make
	make install

## Getting started

In order to test whether PHM is properly installed, run the test script in the *script* directory. It takes some time (should be less than 1 hour) to map regulatory interactions on chromosome 22 using our ATAC-seq data (see *data* directory).

	bash
	export PHMDIR=/path/to/phmdir/
	cd $PHMDIR
	sh script/test.sh

End of the script, you will find the two directories *Stage1* and *Stage2* in the *data* directory, which contain output files of the hierarchical model and pairwise hierarchical model, respectively. The detailed workflow is as follows.

## Workflow

The mapping procedure is split into two stages to reduce the computational complexity of parameter estimation.

![workflow](https://github.com/natsuhiko/Images/blob/master/workflow.png)

### 1st stage

The hierarcical model (HM) takes Bayes factors (BFs) of QTL associations as an input data. We first describe how to compute BFs from normalised read counts and genotype data in a VCF file. The following script *bayeslm1.sh* computes BFs across all variants in the cis-window for each peak.

	sh $PHMDIR/script/bayeslm1.sh 1 2219 \
		/path/to/your/normalized_count.bin \
		/path/to/your/VCF.tbx.gz \
		/path/to/your/peak_bed.tbx.gz \
		/path/to/your/output_file_bayeslm1.gz

The first 2 arguments (1 and 2219) suggests it computes BFs for all 2,219 peaks on chromosome 22 as a single batch. You can split peaks into equaly sized bins for parallelise the job. For example,

	sh $PHMDIR/script/bayeslm1.sh 4 100 ...

to compute the 4th bin with 100 peaks for each bin. The normalized count data can be any quantity at each peak properly normalized so that you can assume normality (e.g., log FPKM, log TPM, quantile normalized value and so on). The file must be formaed as a binary double array. We also provide a simple R script to convert a raw read count table into the specific format we use (see below). The VCF file and the peak annotation files are needed to be tabix indexed. The last argument is the file name of output which is always gzipped. Note that the order of arguments matters for the script.

Once you obtained the BFs, you can feed the output file into the hierarnical model to estimate the variant-level and peak-level prior probabilities: 

	$PHMDIR/bin/hm \
		-i /path/to/your/output_file_bayeslm1.gz \
		-c I,S,S,S,S,S,S,S,C2,C3,S,S,B \
		-r 6186092 -f 2219 \
		-p -o /path/to/your/output_directory_hm

The **-c** option specifies which column of *output_file_bayeslm1.gz* is used as covariates. 

|:---|:-----------|
| I  | Peak ID |
| Cn | Categorical variable with *n* levels |
| Nm | Numerical variable with *m* spline basis |
|    | If *m=0* then it is used as linear |
| B  | Bayes factors |
| S  | Skipped and unused in hm |

Note that the BED file of the peak annotation has to be tabix indexed. The regional Bayes factors for PHM also calculated with the parameter esitimate from the hierarchical model.

	bayeslm -g /path/to/your/VCF/file.gz -w 1000000 -j 10 -f /path/to/your/peak.bed.gz -p1 variant_level_prior.gz -p2 peak_level_prior.gz

The function provides the regional bayes factors for all peak k (>j) in the cis-window.

### 2nd stage

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

## Warnings

Current hm only accepts the ID field to be sequencial and starting from 1. It is not allowed to use, for example, 2,3,4,...(starting from 2) or 1,2,4,...(a gap at 3). 
