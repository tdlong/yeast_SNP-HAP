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
```

This approach is specific to building SNP tables from poolseq samples (or diploid or haploid individuals) obtained from a synthetic population with known and sequenced "founders".  Since the founders are highly chracterized, SNPs can be divided into 3 flavors: 1) all SNPs indentified in the founders (including ones that are not genotyped reliably), 2) well-behaved or what we call Granny SNPs (these are SNPs that can be genotyped very reliably, or what you grandmother would classifly as a SNP based on how you defined a SNP to her, think aggressive filtering), and candidate new mutations (SNPs never seen in the founders, even as bad looking SNPs).  Viewing the world this way is justified for synthetic populations and avoids the GATK shit-show associated with "adding one sample" when you have a large number of samples. We are also finding GATK's newer algorithms (v4) increasingly problematic for poolseq data where the sample is not diploid. So although GATK is the industry standard we only use it to identify and filter SNPs in the founders (what it is designed to do), but then we avoid its use for calling SNPs in samples pulled from the population.  This pipeline assumes you show up to play with the two lists of SNPs obtained from the founders using GATK (but we do not describe that pipeline).

```R

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
library(limSolve)

# https://rdrr.io/cran/limSolve/man/lsei.html
# https://cran.r-project.org/web/packages/limSolve/vignettes/limSolve.pdf

#In matrix notation, linear inverse problems are defined as:
#A·x ~= b  (1)
#E·x = f   (2)
#G·x >= h  (3)
#There are three sets of linear equations:  equalities that have to be met as closely as possible (1), equalities that have to be met exactly (2) and inequalities (3).

# Typically:
# (1) is a linear regression model (A = design matrix {n,d}, b = observations {n}, x = predictors {d})
# (2) is an exact constraint (e.g., the coefficients sum to 1), then E is {ec,d}, f is {ec}, where ec is the number of exact constraints.
# (3) is an inequality (e.g., all coefficients are positive), then G is {ic,d}, h is {ic}, where ic is the number of inequality constraints

# in the example below the xi's sum to one (so E is row matrix of d ones).
# and
# each xi is >= 0 (so G is an Identity matrix of rank d AND H is a column vector of zeros)
# the constraint on the sum means I do not have specify each xi is <= 1 (equivalent to saying each -xi >= -1)

d = ncol(predictors)
A = predictors
B = Y
E = t(matrix(rep(1,d)))
F = 1
G = diag(rep(1,d))
H = matrix(rep(0,d))
Wa = weights          # of course I was also doing a weighted regression ... why not

out = lsei(A=A,B=B,E=E,F=F,G=G,H=H,Wa=Wa,verbose=TRUE)
