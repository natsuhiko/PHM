git add bayeslm.c  getLogBF.c  hm.c phm.c  util_bayeslm.c  util.c
git add bayeslm.h  getLogBF.h  hm.h makeUsage.sh  phm.h  usage.h  util_bayeslm.h  util.h
git add usage_bayeslm.txt  usage.h  usage_hm.txt  usage_phm.txt ../README.md
git add ../script/bayeslm1.sh ../script/bayeslm2.sh ../script/test.sh ../script/test_coloc.sh ../script/bayeslm_eqtl.sh ../script/bayeslm_coloc.sh
git add ../data/variant_level0.bin ../data/genes.chr22.bed.gz ../data/log.fpkm.geu.372.chr22.bin ../data/chr22.geu.372.vcf.gz ../data/chr22.geu.372.vcf.gz.tbi
git commit -m $1
git push origin master
