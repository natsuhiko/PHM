#!/bin/bash

grep -v "^#" usage_phm.txt | sed 's/\\$/\\\\/g' | awk -F '\n' 'BEGIN{print "#include<stdio.h>\nvoid usage(){"}{print "\tfprintf(stderr, \""$1"\\n\");"};END{print "}"}' > usage_phm.c
grep -v "^#" usage_bayeslm.txt | sed 's/\\$/\\\\/g' | awk -F '\n' 'BEGIN{print "#include<stdio.h>\nvoid usage(){"}{print "\tfprintf(stderr, \""$1"\\n\");"};END{print "}"}' > usage_bayeslm.c
