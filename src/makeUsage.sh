#!/bin/bash

grep -v "^#" usage_phm.txt | sed 's/\\$/\\\\/g' | awk -F '\n' 'BEGIN{print "#include<stdio.h>\n\nvoid usage_phm(){"}{print "\tfprintf(stderr, \""$1"\\n\");"};END{print "}"}' > usage.c
grep -v "^#" usage_hm.txt | sed 's/\\$/\\\\/g' | awk -F '\n' 'BEGIN{print "\nvoid usage_hm(){"}{print "\tfprintf(stderr, \""$1"\\n\");"};END{print "}"}' >> usage.c
grep -v "^#" usage_bayeslm.txt | sed 's/\\$/\\\\/g' | awk -F '\n' 'BEGIN{print "\nvoid usage_bayeslm(){"}{print "\tfprintf(stderr, \""$1"\\n\");"};END{print "}"}' >> usage.c
