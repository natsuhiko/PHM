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

The first 2 arguments (1 and 2219) suggest it computes BFs for all 2,219 peaks on chromosome 22 as a single batch. You can split peaks into equaly sized bins for parallelise the job. For example,

	sh $PHMDIR/script/bayeslm1.sh 4 100 ...

to compute the 4th bin with 100 peaks for each bin. The normalized count data can be any quantity at each peak properly normalized so that you can assume normality (e.g., log FPKM, log TPM, quantile normalized value and so on). The file must be formaed as a binary double array. We also provide a simple R script to convert a raw read count table into the specific format we use (see below). The VCF file and the peak annotation files are needed to be tabix indexed. The last argument is the file name of output which is always gzipped. Note that the order of arguments matters for the script.

Once you obtained the BFs, you can feed the output file into the hierarnical model to estimate the variant-level and peak-level prior probabilities: 

	$PHMDIR/bin/hm \
		-i /path/to/your/output_file_bayeslm1.gz \
		-c I,S,S,S,S,S,S,S,C2,C3,S,S,B \
		-r 6186092 -f 2219 \
		-p -o /path/to/your/output_directory_hm

The **-c** option specifies which column of *output_file_bayeslm1.gz* is used as input for hm. In this example, the first column is the ID of peaks, the 9th column is variant type as a categorical variable (0=SNP; 1=INDEL; 2=CNV), the 10th column is location of variants as a categorical varialbe (0=outsize peak; 1=inside the focal peak; 2=inside a flanking peak) and the 13th column is BFs. Other columns are ignored in the model fitting. The table below illustrates each column type and its description. <span style="color:red">__Note that, for pairwise hierarhical model, the above column specification is only vaiable for current implementation.__</span>

| Column Type | Description | 
|:----:|:-----------------------------------------|
| I    | Peak ID                                  |
| C*n* | Categorical variable with *n* levels. If *n* is smaller than the actual number of levels, levels >*n* is treated as the level *n*. |
| N*m* | Numerical variable with *m* spline bases. The variable must be scaled in [0,1]. If *m*=0, then the variable is used as a linear predictor (not necessarily scaled in this case). |
| B    | Bayes factors                            |
| S    | Skipped and unused in hm                 |

### 2nd stage

Pairwise hierarchical model takes regional Bayes factors (RBFs) which are calculated by using the first stage parameter estimate. The following script *bayeslm2.sh* computes RBFs for any peak pairs in the cis-window:

	sh $PHMDIR/script/bayeslm2.sh 1 2219 \
		/path/to/your/normalized_count.bin \
		/path/to/your/VCF.tbx.gz \
		/path/to/your/peak_bed.tbx.gz \
		/path/to/your/output_directory_hm/variant_level.bin \
		/path/to/your/output_directory_hm/Pi1.bin \
		/path/to/your/output_file_bayeslm2.gz

Again, the first 2 arguments suggest the script computes RBFs for all peaks as a single batch. You can parallelise the script as same as before. The 6th and 7th arguments specify the output parameter estimate from the hierarchical model (hm) in the 1st stage. The 8th argument is the output RBFs in gzipped format. 

After computing RBFs, you can fit the pairwise hierarhical model (*phm*) to estimate the peak-pair-level prior probability:

	$PHMDIR/bin/phm \
		-i /path/to/your/output_file_bayeslm1.gz \
		-c J,K,N4,B10 \
		-r 85323 -f 2219 \
		-o /path/to/your/output_directory_phm

The **-c** option specifies which column of *output_file_bayeslm2.gz* is used as input for phm. In this example, the first 2 columns are the peak ID for 5' and 3' ends,
the 3rd column is the peak distance between paired peaks and the 4th column is the RBFs for the 10 different interaction sub-hypotheses (in the space separated format). The table below illustrates each column type and its description. 

| Column Type | Description |
|:----:|:-----------------------------------------|
| J    | Peak ID at 5' end                                |
| K    | Peak ID at 3' end                                |
| C*n* | Categorical variable with *n* levels. If *n* is smaller than the actual number of levels, levels >*n* is treated as the level *n*. |
| N*m* | Numerical variable with *m* spline bases. The variable must be scaled in [0,1]. If *m*=0, then the variable is used as a linear predictor (not necessarily scaled in this case). |
| B10  | Regional Bayes factors for the 10 interaction sub-hypotheses                |
| S    | Skipped and unused in phm                 |

As outputs, you will obtain the posterior probability of each interaction hypothesis (pp.gz) and the probability of master regulator (pmr.gz).
The following shows the column descriptions:

| Column No. | Description (pp.gz) |
|:----:|:-----------------------------------------|
| 1    | Peak ID at 5' end                              |
| 2    | Peak ID at 3' end                              |
| 3    | Two peaks are not QTL (null) |
| 4    | 5' Peak is QTL (single) |
| 5    | 3' Peak is QTL (single) |
| 6    | Two peaks are indelendent QTLs (linkage) |
| 7    | Two peaks are in pleiotoropic (pleiotropy) |
| 8    | 5' peak regulates 3' peak (causality) |
| 9    | 3' peak regulates 5' peak (causality) |

| Column No. | Description (pmr.gz) |
|:----:|:-----------------------------------------|
| 1    | Probability of master regulator (PMR) given that the peak is a QTL |
| 2    | Expected number of downstream peaks                             |
| 3    | Expected number of upstream peaks |
| 4    | Maximum a posteriori parent peak ID |
| 5    | Posterior probability of causality for the MAP parent |

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
