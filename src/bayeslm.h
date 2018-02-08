#include <config.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include "htslib/tbx.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/regidx.h"

#include <zlib.h>
#include <math.h>
#include "getLogBF.h"

int exp_gt_gtdsgl;
int verbose_loadVCF;


#define FORMAT_GT 1
#define FORMAT_GL 2
#define FORMAT_AP 3
#define FORMAT_GP 4
#define FORMAT_AS 5
#define FORMAT_RD 6
#define FORMAT_BF 7
#define FORMAT_DS 8
#define FORMAT_OTHER 0


#define VT_SNP 0
#define VT_INDEL 1
#define VT_SV 2
#define VT_OTHER 3

typedef struct{
        int VT;
        double RSQ;
        double AF;
        double CER;
}VCF_info;


#define MODE_N 0
#define MODE_P 1
#define MODE_C 2
#define MODE_G 3
#define MODE_PP 4


