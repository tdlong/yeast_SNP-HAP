# yeast_SNP-HAP
pipeline for calling SNPs and Haplotypes, associated with resource paper
```bash
########################
# programs needed
########################

# vcflib and vt are not as common as the others

vcflib       # https://github.com/vcflib/vcflib
vt/v0.57721  # http://genome.sph.umich.edu/wiki/vt
R/3.4.1
bcftools/1.3
enthought_python/7.3.2
bwa/0.7.8
samtools/1.9
freebayes/0.9.21 
vcftools/0.1.15 
#  bedtools/2.25.0 creates a conflict with vt & vcfallelicprimitives
bedtools/2.15.0

########################
#  make directories
########################

mkdir helpfiles
mkdir scripts
mkdir May1
mkdir raw
cd May1
mkdir SNPs
mkdir haplos
mkdir bam
mkdir mut
cd haplos
mkdir T3

#######################
# helpfiles and scripts needed to run
#######################

helperfiles/
├── allFounders.SNPlist
├── barcodes.txt
├── chromes.longfirst.txt
├── files.May1.txt
├── founder.file.sexual.June11.txt
├── Granny.in.allFound.txt
└── sexual.names.txt

scripts/
├── haplotype.limSolve.general.sh
├── haplotyper.limSolve.code.R
├── haplotyper.limSolve.R
├── haplotyper.merge.sh
├── magic_count.pl
├── make.muttable.general.pl
├── make.SNPtable.general.pl
├── process_dipx.sh
└── rename.reads.general.py


########################
## process fastq files ...
########################

qsub scripts/process_dipx.sh

########################
## exclude failed samples...
########################

cd May1/bam
ls -l *.bam | awk '{if ($5 > 10000) print $9}' | sed 's/.bam//'  >> .../helperfiles/files.May1.txt
# I excluded these as "Nextera fails", note the file sizes of the bams!!!
ls -l *.bam | awk '{if ($5 <= 10000) print}' 

########################
##  SNPs
########################

input=".../helperfiles/files.May1.txt"
perl scripts/make.SNPtable.general.pl "May1/SNPs/" ".known.adl" <$input >SNPtable.May1.txt
# sort of first two fields...
cat SNPtable.May1.txt | head -n 1 > SNPtable.May1.sort.txt
cat SNPtable.May1.txt | tail -n +2 | sort -k 1,1 -k 2,2n >> SNPtable.May1.sort.txt
#check
cat SNPtable.May1.sort.txt | cut -f 1,2,3,5,7,9,11,13,15,17,19,43 | more 
rm SNPtable.May1.txt

##########
##  new mutations ...
##########

perl scripts/make.muttable.general.pl "May1/mut/" ".newmut.vcf" <$input >temp.txt
head -n 1 temp.txt > MUT.dipx.May1.txt
tail -n +2 temp.txt | sort -k1,1 -k2,2n -k3,3 >> MUT.dipx.May1.txt
rm temp.txt

##########
##  make list of samples and chromosomes to allow parallel processing ...
##########
module load R/3.5.3 
R
samples=scan(file=".../helperfiles/files.May1.txt",what=character())
chromes=scan(file=".../helperfiles/chromes.longfirst.txt",what=character())
for (i in chromes){
	for (j in samples){
		cat(file="chr_pool_May1.txt",i,"\t",j,"\n",append=TRUE)
		}
	}

wc -l chr_pool_May1.txt
# 7242

##########
## run haplotyper ...
##########

input="chr_pool_May1.txt"
qsub -t 1-7242 scripts/haplotype.limSolve.general.sh $input May1/haplos/T3 SNPtable.May1.sort.txt helperfiles/founder.file.sexual.June11.txt

# wait until done (check all samples were processed -> ls | wc -l)
# sh haplotyper.merge.sh dir get_header outfile
sh scripts/haplotyper.merge.sh May1/haplos/T3 BAS02_chrI_hap_freq.txt haps.dipx.June11.F15AW5050.limSolve.txt.gz

```
