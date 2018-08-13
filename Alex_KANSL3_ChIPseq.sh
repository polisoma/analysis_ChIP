mkdir ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/raw

# from my Desktop:
# scp -r ChipseqARD01-* daria@tycho.sund.root.ku.dk:~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/raw/.

# pulling lanes together
 
for f in `ls raw/Chipseq*/*_R1_001.fastq.gz`; do   
   name=$(echo -en "${f%%.*}" | awk '{split($1,a,"_S"); print a[1]}')
   label=$(basename $name)
   echo -en "$label\n"
   zcat $f >> ~/temp_NOBACKUP/$label.merged.fastq
done

# number of reads: 

for name in T-KANSL3 K-KANSL3 T-INPUT K-INPUT; do
   echo -en "$name\n"
   cat ~/temp_NOBACKUP/$name.merged.fastq | wc -l
done

# 122006360/4=30501590
# 162299576/4=40574894
# 211203700/4=52800925
# 192335536/4=48083884

screen -R mapping
# do not have gz files - cat on files 
# will map to hg19 as all Fantom data is for hg19 + motifs
hg19=/k/genomes/hg19/index/bowtie_canonical/hg19 # somehow index in my directory was not working
for name in T-KANSL3 K-KANSL3 T-INPUT K-INPUT; do
  echo -en $name
  cat ~/temp_NOBACKUP/$name.merged.fastq | bowtie -p 5 -q -m 1 -v 3 --sam --best --strata --quiet $hg19 - > ~/temp_NOBACKUP/$name.sam
done

for sample in T-KANSL3 K-KANSL3 T-INPUT K-INPUT; do
  echo $sample
  samtools view -Sb ~/temp_NOBACKUP/${sample}.sam > ~/temp_NOBACKUP/${sample}_nonSorted.bam
  samtools sort ~/temp_NOBACKUP/${sample}_nonSorted.bam > ~/temp_NOBACKUP/${sample}.bam
  samtools index ~/temp_NOBACKUP/${sample}.bam
done

# after everything is done
for sample in T-KANSL3 K-KANSL3 T-INPUT K-INPUT; do
rm ~/temp_NOBACKUP/${sample}.sam ~/temp_NOBACKUP/${sample}_nonSorted.bam
done

for sample in T-KANSL3 K-KANSL3 T-INPUT K-INPUT; do
mv ~/temp_NOBACKUP/${sample}.bam  ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/bam/.
mv ~/temp_NOBACKUP/${sample}.bam.bai  ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/bam/.
done 

# making bw files
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017
for sample in T-KANSL3 K-KANSL3 T-INPUT K-INPUT; do
   echo $sample
  EXTEND=140
  CHROMSIZE=/home/daria/data/data_genomes/hg19/chrom.sizes
  # Number of reads
  librarySize=$(samtools idxstats bam/${sample}.bam | awk '{total+=$3} END{print total}')
  echo $librarySize
  # Create density file: extend reads, calculate read density at each position and normalize the library size to 1 million reads
  bamToBed -i bam/${sample}.bam | awk -vCHROM=$CHROMSIZE -vEXTEND=$EXTEND -vOFS='\t' 'BEGIN{while((getline<CHROM)>0){chromSize[$1]=$2}}{chrom=$1;start=$2;end=$3;strand=$6;if(strand=="+"){end=start+EXTEND;if(end>chromSize[chrom]){end=chromSize[chrom]}};if(strand=="-"){start=end-EXTEND;if(start<1){start=1}};print chrom,start,end}' | sort -k1,1 -k2,2n | genomeCoverageBed -i stdin -g $CHROMSIZE -d | awk -vOFS='\t' -vSIZE=$librarySize '{print $1,$2,$2+1,$3*1000000/SIZE}' | gzip > ~/temp_NOBACKUP/${sample}.density.gz	
  # Create WIG file
  gunzip -c ~/temp_NOBACKUP/${sample}.density.gz | awk -vOFS='\t' '($4!=0){if(!chrom[$1]){print "variableStep chrom="$1;chrom[$1]=1};print $2,$4}' | gzip > ~/temp_NOBACKUP/${sample}.wig.gz
  # Create BigWig file, takes quite some time and memory-consuming
  wigToBigWig ~/temp_NOBACKUP/${sample}.wig.gz $CHROMSIZE bw/${sample}.bw
  # Remove intermediate file
  #rm ${sample}.wig.gz ${sample}.density.gz
done

######################
# putting files in respective folders
mkdir ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/bam
mkdir ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/bw


###################### calling peaks
# macs
# following this: https://github.com/taoliu/MACS/wiki/Advanced:-Call-peaks-using-MACS2-subcommands
mkdir ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/macs2
screen -R macs
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/

for pair in T-KANSL3_T-INPUT K-KANSL3_K-INPUT; do
   exp=$(echo $pair | awk '{split($1,a,"_");print a[1]}')
   input=$(echo $pair | awk '{split($1,a,"_");print a[2]}')
   echo -en "$exp\t$input"
   macs2 callpeak -t bam/$exp.bam -c bam/$input.bam -f BAM -n ${exp}_${input} --outdir macs2 --call-summits -g hs --keep-dup 1
done

# filtering dup and making bed files
for sample in T-INPUT K-INPUT; do
   macs2 filterdup -i bam/${sample}.bam --keep-dup=1 -o macs2/${sample}_filterdup.bed
done&

# random extsize
# making pileup in macs2
# shoudl use 140 extension
#macs2 pileup -i macs2/K-KANSL3_filterdup.bed -o macs2/K-KANSL3_filterdup.pileup.bdg --extsize 130 &
#macs2 pileup -i macs2/T-KANSL3_filterdup.bed -o macs2/T-KANSL3_filterdup.pileup.bdg --extsize 130 &

# bedGraphtoBigWig
# claims that chrM is bigger
# removed chrM
CHROMSIZE=chrom.sizes # problems with chromoc=somal sizes
bedGraphToBigWig macs2/K-KANSL3_filterdup.pileup.bdg $CHROMSIZE bw/K-KANSL3_filterdup.pileup.bw
bedGraphToBigWig macs2/T-KANSL3_filterdup.pileup.bdg $CHROMSIZE bw/T-KANSL3_filterdup.pileup.bw

#

# adding the first line to narrowPeak
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017
(echo -en "track type=narrowPeak visibility=3 db=hg19 name=\"K-KANSL3_K-INPUT_narrow_peak\" description=\"K-KANSL3_K-INPUT_narrow_peak\"\n" ; cat macs2/K-KANSL3_K-INPUT_peaks.narrowPeak ) > macs2/K-KANSL3_K-INPUT_peaks.browser.narrowPeak
(echo -en "track type=narrowPeak visibility=3 db=hg19 name=\"T-KANSL3_T-INPUT_narrow_peak\" description=\"T-KANSL3_T-INPUT_narrow_peak\"\n" ; cat macs2/T-KANSL3_T-INPUT_peaks.narrowPeak ) > macs2/T-KANSL3_T-INPUT_peaks.browser.narrowPeak

# loading peaks and bw to the browser
cd ~/web/projects/Alex
ln -s /NextGenSeqData/project-data/daria/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/macs2/K-KANSL3_K-INPUT_peaks.browser.narrowPeak K-KANSL3_K-INPUT_peaks.browser.narrowPeak
ln -s /NextGenSeqData/project-data/daria/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/macs2/T-KANSL3_T-INPUT_peaks.browser.narrowPeak T-KANSL3_T-INPUT_peaks.browser.narrowPeak
ln -s /NextGenSeqData/project-data/daria/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/bw/*.bw .

# links
https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/K-KANSL3_K-INPUT_peaks.browser.narrowPeak
https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/T-KANSL3_T-INPUT_peaks.browser.narrowPeak

https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/diff_c1_vs_c2_c3.0_cond2.browser.bed

# tracks
track type=bigWig name="K-KANSL3.bw" description="K-KANSL3" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/K-KANSL3.bw
track type=bigWig name="T-KANSL3.bw" description="T-KANSL3" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/T-KANSL3.bw
track type=bigWig name="T-INPUT.bw" description="T-INPUT" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/T-INPUT.bw
track type=bigWig name="K-INPUT.bw" description="K-INPUT" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/K-INPUT.bw

##### differential peaks with macs2: https://github.com/taoliu/MACS/wiki/Call-differential-binding-events
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017
mkdir -p ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/macs2/differential

macs2 callpeak -B -t bam/K-KANSL3.bam -c bam/K-INPUT.bam -n K-ChIP --outdir macs2/differential --nomodel --extsize 200 &
macs2 callpeak -B -t bam/T-KANSL3.bam -c bam/T-INPUT.bam -n T-ChIP --outdir macs2/differential --nomodel --extsize 200 &

# differential peak calling
# need to know number of tags after filtering
# total tags in treatment after filtering: 32884992; in input: 39818748
# total tags in control after filtering: in control 43796969

tags1=39818748 # for K
tags2=43796969 # for T
folder=macs2/differential
# -g and -l parameters left unchanged
macs2 bdgdiff --t1 $folder/K-ChIP_treat_pileup.bdg --c1 $folder/K-ChIP_control_lambda.bdg --t2 $folder/T-ChIP_treat_pileup.bdg\
   --c2 $folder/T-ChIP_control_lambda.bdg --d1 $tags1 --d2 $tags2 -g 60 -l 120 --o-prefix $folder/diff_c1_vs_c2

######### calling peaks with macs2, no model
# see above for differential peaks --> 68 peaks  

######################################### peakzilla
mkdir -p ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/peakzilla
peakzilla=~/utils/peakzilla/peakzilla-master/peakzilla.py

# bed files for peakzilla were generated with macs2
#for sample in T-INPUT K-INPUT; do
#   macs2 filterdup -i bam/${sample}.bam --keep-dup=1 -o macs2/${sample}_filterdup.bed
#done&

$peakzilla macs2/T-KANSL3_filterdup.bed macs2/T-INPUT_filterdup.bed > peakzilla/T-KANSL3_peakzilla.tsv
$peakzilla macs2/K-KANSL3_filterdup.bed macs2/K-INPUT_filterdup.bed > peakzilla/K-KANSL3_peakzilla.tsv

$peakzilla macs2/T-KANSL3_filterdup.bed > peakzilla/T-KANSL3_peakzilla.noinput.tsv
$peakzilla macs2/K-KANSL3_filterdup.bed > peakzilla/K-KANSL3_peakzilla.noinput.tsv

###################### genomic distribution
#wget ftp://ftp.ensembl.org/pub/release-71/gtf/homo_sapiens/Homo_sapiens.GRCh37.71.gtf.gz
#gunzip -c Homo_sapiens.GRCh37.71.gtf.gz | awk -vOFS='\t' '($2=="protein_coding"&&(($1>=1&&$1<=22)||$1=="Y"||$1=="X")){print "chr"$0}' | sort -k1,1 -k4,4n > Homo_sapiens.GRCh37.71.gtf
#gunzip -c Homo_sapiens.GRCh37.71.gtf.gz | awk -vFS='\t' -vOFS='\t' -vC="/work2/gschub/anais/data/chrom/hg19.chrom.sizes.good" 'BEGIN{while(getline<C){S[$1]=$2}}{if("chr"$1 in S && $2=="protein_coding"){gid=$9;sub(/gene_id "/,"",gid);sub(/";.*/,"",gid);tid=$9;sub(/.*transcript_id "/,"",tid);sub(/";.*/,"",tid);gn=$9;sub(/.*gene_name "/,"",gn);sub(/";.*/,"",gn);print "chr"$1,$4,$5,$7,$3,gid,tid,gn}}' | sort -k1,1 -k2,2n | gzip > hg19_annotation.txt.gz

# annotations are in ~/data/data_genomes/annotations/raw/
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/
fileA=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
fileB=~/data/data_genomes/annotations/raw/hg19_genomic_features.txt

bedtools intersect -wb -a <(awk -vOFS="\t" '{print $1,$2+int(($3-$2)/2),$2+int(($3-$2)/2)}' $fileA | sort -k1,1 -k2,2n) -b <(awk -vOFS="\t" '(NR>1){print $0}' $fileB | sort -k1,1 -k2,2n) | awk '{count[$7]++}END{for (i in count)print i,"\t"count[i]}'
# for the whole genome
awk '{count[$4]++}END{for (i in count)print i,"\t"count[i]}' $fileB

# overlaping with CPG islands
fileA=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
fileB=~/data/data_genomes/annotations/hg19_cgi.txt

bedtools intersect -c -a <(awk -vOFS="\t" '{print $1,$2+int(($3-$2)/2),$2+int(($3-$2)/2)}' $fileA | sort -k1,1 -k2,2n) -b <(awk -vOFS="\t" '(NR>1){print $0}' $fileB | sort -k1,1 -k2,2n) | awk '($4==1){count++}END{print count}'

# overlapping with FANTOM peaks

# preprocessing of FANTOM peaks
#(echo -en "chr\tstart\tend\tstrand\tgene_FANTOM\tglobal_promoter_rank\tCAGE_fresh\tCAGE_revived\tCAGE_thawed\n"
#awk -vOFS="\t" '(NR>1){gsub(/\../,"=",$1) ;a=substr($1,0,length($1)-2); split($1,strand,","); split(a,chr,":"); split(chr[2],pos,"="); #if($4!="NA"){split($2,temp,"@");split(temp[2],g,",");gene=g[1];rank=temp[1]}else{gene="NA";rank="NA"} ; print #chr[1],pos[1],pos[2],strand[2],gene,rank,$7,$8,$9}' hg19.cage_peak_tpm_ann_decoded.osc.txt.gz.extract_THP1.tsv ) > FANTOM_table_THP1.txt

# the same for THP1 but no mean, taking only CAGE_fresh into account
cd /NextGenSeqData/project-data/daria/projects/BRIC_data/Alex/CRISPR_CAGE_files/data
awk -vOFS="\t" '(NR>1){print $0}' FANTOM_table_THP1.txt | sort -k5,5 -k7,7gr | awk '{if($5!="NA"){if($7!=0){if(NR==1 || $5!=gene){n=1;gene=$5;print $0, "pRank"n}else{n++; print $0,"pRank"n}}else{print $0,0}}else{print $0,"NA"}}' > ~/temp_NOBACKUP/FANTOM_table_THP1.rank.txt

#bedintersect somehow is not so good
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/

fileA=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
fileB=~/temp_NOBACKUP/FANTOM_table_THP1.rank.txt
bedtools intersect -loj -a <(awk -vOFS="\t" '{print $1,$2,$3}' $fileA | sort -k1,1 -k2,2n) -b <(awk -vOFS="\t" '(NR>1){print $0}' $fileB | sort -k1,1 -k2,2n) | awk '($4=="."){print $0}' | wc -l 

# closest FANTOM5 peaks --> best solution, better than intersect
# if overlaps several --> will report all
# temp file does not have a header! mistake!!!!
# NO HEADER! 
bedtools closest -a <(awk -vOFS="\t" '{print $1,$2,$3}' $fileA | sort -k1,1 -k2,2n) -b <(awk -vOFS="\t" '(NR>1){print $0}' $fileB | sort -k1,1 -k2,2n) -d | awk '($8=="NA")'

bedtools closest -a <(awk -vOFS="\t" '{print $1,$2,$3}' $fileA | sort -k1,1 -k2,2n) -b <(awk -vOFS="\t" '(NR>1){print $0}' $fileB | sort -k1,1 -k2,2n) -d | awk '{print $8}' | uniq > /tmp/genes.txt

# overlapping CAGE peaks: distribution compare to other CAGE peaks of expressed genes
fileA=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
fileB=~/temp_NOBACKUP/FANTOM_table_THP1.rank.txt

awk '{print $7}' $fileB | median.R -i - -n -q # most of the genes are not expressed
awk '($7>1){print $7}' $fileB | median.R -i - -n -q # median 5.22
awk '($7!=0){print $7}' $fileB | median.R -i - -n -q $ median 1.23

# $10 column - CAGE from fresh THP1
# distance shoudl be <=200 --> consider as an overlap
bedtools closest -a <(awk -vOFS="\t" '{print $1,$2,$3}' $fileA | sort -k1,1 -k2,2n) -b <(awk -vOFS="\t" '{print $0}' $fileB | sort -k1,1 -k2,2n) -d | awk '($14<200){print $10}' | median.R -i - -n -q # median 4.45 -- highly expressed

bedtools closest -a <(awk -vOFS="\t" '{print $1,$2,$3}' $fileA | sort -k1,1 -k2,2n) -b <(awk -vOFS="\t" '{print $0}' $fileB | sort -k1,1 -k2,2n) -d | awk '($14<200){print $10}' > /tmp/overlap.txt

# random CAGE peaks and values in CAGE fresh THP1
awk '($7!=0){print $7}' $fileB  | shuf -n 1000 > /tmp/random.txt

# making plots in R
R
source("~/utils/Rutils/boxplot95.R")

t1<-read.table("~/temp_NOBACKUP/FANTOM_table_THP1.rank.txt")
t2<-read.table("/tmp/overlap.txt")
t3<-read.table("/tmp/random.txt")

par(mfrow=c(1,3)) 
boxplot95(log(t1[t1[,7]!=0,7])/log(2),frame=F,ylim=c(-3,6))
boxplot95(log(t2)/log(2),frame=F,ylim=c(-3,6))
boxplot95(log(t3)/log(2),frame=F,ylim=c(-3,6))

q()

# confirmation CAGE observation with RNA-seq data
# taking genes where peaks are <200 bp from CAGE

fileA=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
fileB=~/temp_NOBACKUP/FANTOM_table_THP1.rank.txt
bedtools closest -a <(awk -vOFS="\t" '{print $1,$2,$3}' $fileA | sort -k1,1 -k2,2n) -b <(awk -vOFS="\t" '{print $0}' $fileB | sort -k1,1 -k2,2n) -d | awk '($14<200){print $8}' | uniq > /tmp/genes.txt

awk -vOFS="\t" -vFile="/tmp/genes.txt" 'BEGIN{while((getline<File)>0){gene[$1]=$1}}{if($1 in gene){print $1,($15+$16+$17)/3}}' ../../RNA-seq/data_May2017/RPKMs/all_counts_rpkm.featureCount.final.txt | cut -f2 | median.R -i - -q -n # 0  4.4182  median=8.48329 16.65335  530.967

awk '(($15+$16+$17)/3!=0){print ($15+$16+$17)/3}' ../../RNA-seq/data_May2017/RPKMs/all_counts_rpkm.featureCount.final.txt | median.R -i - -q -n #median=0.511658

# should also take closest TSS instead of relying on CAGE data assignment

R
source("~/utils/Rutils/boxplot95.R")

t1<-read.table("/tmp/rna2.txt")
t2<-read.table("/tmp/rna1.txt")

par(mfrow=c(1,3)) 
boxplot95(log(t1[,1])/log(2),frame=F,ylim=c(-7,6))
boxplot95(log(t2[,1])/log(2),frame=F,ylim=c(-7,6))

q()

# overlap of genes with RNA-seq
# files are in ../../RNA-seq/data_May2017/DESeq2/

awk -vOFS="\t" -vFile="/tmp/genes.txt" 'BEGIN{while((getline<File)>0){gene[$1]=$1}}{if($2 in gene){print $2,$3,$4,$5,$6,$7,$8}}' ../../RNA-seq/data_May2017/DESeq2/all_experiemnts.deseq2.merged.txt | awk '($3<=0.01 || $5<=0.01 || $7<=0.01)' | head

awk -vOFS="\t" -vFile="/tmp/genes.txt" 'BEGIN{while((getline<File)>0){gene[$1]=$1}}{if($2 in gene){print $2,$3,$4,$5,$6,$7,$8}}' ../../RNA-seq/data_May2017/DESeq2/all_experiemnts.deseq2.merged.txt | awk '($3<=0.01 || $5<=0.01 || $7<=0.01){print $1}' > /tmp/overlap.genes.txt

awk -vOFS="\t" '($3<=0.01 || $5<=0.01 || $7<=0.01){print $2,$3,$4,$5,$6,$7,$8}' ../../RNA-seq/data_May2017/DESeq2/all_experiemnts.deseq2.merged.txt > /tmp/diffregulated.txt

# the changes were not so dramatic --> will look at the distribution
awk -vOFS="\t" '{print $2,$3,$4,$5,$6,$7,$8}' ../../RNA-seq/data_May2017/DESeq2/all_experiemnts.deseq2.merged.txt | awk '($7<=0.01){print $6}' | median.R -i - -q -n

awk -vOFS="\t" -vFile="/tmp/genes.txt" 'BEGIN{while((getline<File)>0){gene[$1]=$1}}{if($2 in gene){print $2,$3,$4,$5,$6,$7,$8}}' ../../RNA-seq/data_May2017/DESeq2/all_experiemnts.deseq2.merged.txt | awk '($7<=0.01){print $6}' | median.R -i - -q -n

#R
#library(ggplot2)
#t<-read.table("/tmp/diffregulated.txt")

#par(mfrow = c(1,3))
#for (n in c(2,4,6)){
#genes <- scan("/tmp/overlap.genes.txt",character(), quote = "")
#boxplot(t[,n],ylim=c(-4,1),frame=FALSE)
#for (i in 1:length(genes)){
#y<-c(t[grep(paste0("^",genes[i],"$"),t[,1]),n])
#points<-c(rep(1,t[grep(paste0("^",genes[i],"$"),t[,1]),n]))
#lines(jitter(points),y,type="p")
#text(t[grep(paste0("^",genes[i],"$"),t[,1]),n], labels=genes[i], cex= 0.7, pos=4,xpd=NA)
#}
#}

#points<-c(rep(1,length(log(dist.Pc.min[,2])/log(10))))
#lines(jitter(points),y,type="p")

###################### De novo peak analysis
# de novo against all CpG islands with DREME

cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/
mkdir -p dreme
fileA=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
g=/home/daria/data/data_genomes/hg19/hg19.fa

# see different parameters in help: bedtools getfasta -h 
bedtools getfasta -fi $g -s -bed <(awk '{print $1"\t"$2"\t"$3"\t"NR"\t""+"}' $fileA) -fo dreme/KANSL3_only_peaks.fa -name
# negative regions: some CAGE peaks of expressed genes (CAGEfresh>1)
fileA=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
fileB=~/temp_NOBACKUP/FANTOM_table_THP1.rank.txt
bedtools intersect -v -a <(awk -vOFS="\t" '($7>=1){print $0}' $fileB | sort -k1,1 -k2,2n) -b <(awk -vOFS="\t" '{print $1,$2,$3}' $fileA | sort -k1,1 -k2,2n)  |  awk '{print $1,$2,$3}' | shuf -n 1000 > dreme/negative.txt
# CAGE peaks are tiny, so will add +/- 100 bp
WIN=100
bedtools getfasta -fi $g -s -bed <(awk -vWIN=$WIN '{print $1"\t"$2-WIN"\t"$3+WIN"\t"NR"\t""+"}' dreme/negative.txt) -fo dreme/negative_win$WIN.fa -name

folder=dreme/KANSL3_v1; mkdir -p $folder
dreme -png -oc $folder -p dreme/KANSL3_only_peaks.fa -n dreme/negative_win100.fa -e 0.05


######################## data from Asifa #################
track type=bigWig name="KANSL3_rep2_mES.bw" description="KANSL3 rep2 mES" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/published_data/GSM1251952_NSL3B_mlES_To1x.bw 
track type=bigWig name="KANSL3_rep1_mES.bw" description="KANSL3 rep1 mES" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/published_data/GSM1251951_NSL3A_mlES_To1x.bw
track type=bigWig name="KANSL3_rep1_NPC.bw" description="KANSL3 rep1 NPC" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/published_data/GSM1251959_NSL3A_mlNPC_To1x.bw
track type=bigWig name="KANSL3_rep2_NPC.bw" description="KANSL3 rep2 NPC" bigDataUrl=https://bricweb.sund.ku.dk/bric-data/daria/published_data/GSM1251960_NSL3B_mlNPC_To1x.bw
 
####### getting fastq data
cd ~/data/published_data/mES/Asifa_KANSL3_mES
mkdir -p fastq

#SRR1258425 KANSL3_mES_rep1
#SRR1258426 KANSL3_mES_rep2
#SRR1258415	input_mES_rep1
#SRR1258416	input_mES_rep2
#SRR1258427 input_NPC_rep1
#SRR1258428 input_NPC_rep2
#SRR1258433 KANSL3_NPC_rep1
#SRR1258434 KANSL3_NPC_rep2

# regular prefetch did nto work
/home/sam/software/sratoolkit.2.8.2-1-ubuntu64/bin/prefetch.2.8.2 -v SRR1258425 SRR1258426 SRR1258415 SRR1258416 SRR1258427 SRR1258428 SRR1258433 SRR1258434

screen -R mapping2
fastqdump=/home/sam/software/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump.2.8.2
cd ~/data/published_data/mES/Asifa_KANSL3_mES
for SRR in SRR1258425 SRR1258426 SRR1258415 SRR1258416 SRR1258427 SRR1258428 SRR1258433 SRR1258434 ; do
   echo -en "$SRR\n"
$fastqdump --outdir fastq --split-files /home/daria/ncbi/public/sra/${SRR}.sra
done

# mapping bowtie PE data
# --strata  do not apply in paired-end mode

cd ~/data/published_data/mES/Asifa_KANSL3_mES
mm9=/k/genomes/mm9/index/bowtie_canonical/mm9
for SRR in SRR1258425 SRR1258426 SRR1258415 SRR1258416 SRR1258427 SRR1258428 SRR1258433 SRR1258434 ; do
   echo -en "$SRR\n"
   bowtie -p 5 -q -m 1 -v 3 --sam --quiet $mm9 -1 fastq/${SRR}_1.fastq -2 fastq/${SRR}_2.fastq > ~/temp_NOBACKUP/${SRR}.sam
done

# had to move the files to clean Tycho
# new location: ~/data/temp_NOBACKUP/

screen -R mapping  #31561.mapping2
cd ~/data/published_data/mES/Asifa_KANSL3_mES
for sample in SRR1258425 SRR1258426 SRR1258415 SRR1258416 SRR1258427 SRR1258428 SRR1258433 SRR1258434; do
  echo $sample
  samtools view -Sb ~/data/temp_NOBACKUP/${sample}.sam > ~/data/temp_NOBACKUP/${sample}_nonSorted.bam
  samtools sort ~/data/temp_NOBACKUP/${sample}_nonSorted.bam > ~/data/temp_NOBACKUP/${sample}.bam
  samtools index ~/data/temp_NOBACKUP/${sample}.bam
done

# after everything is done: move and remove files
mkdir -p ~/data/published_data/mES/Asifa_KANSL3_mES/bam 
folder=~/data/published_data/mES/Asifa_KANSL3_mES/bam 

for sample in SRR1258425 SRR1258426 SRR1258415 SRR1258416 SRR1258427 SRR1258428 SRR1258433 SRR1258434; do
mv ~/data/temp_NOBACKUP/${sample}.bam  $folder/.
mv ~/data/temp_NOBACKUP/${sample}.bam.bai  $folder/.
done 

for sample in SRR1258425 SRR1258426 SRR1258415 SRR1258416 SRR1258427 SRR1258428 SRR1258433 SRR1258434; do
rm ~/data/temp_NOBACKUP/${sample}.sam ~/data/temp_NOBACKUP/${sample}_nonSorted.bam
done

######## peak calling with macs2 

# giving proper names to files and making soft links
mkdir -p bam/links
for sample in SRR1258425 SRR1258426 SRR1258415 SRR1258416 SRR1258427 SRR1258428 SRR1258433 SRR1258434; do
  name=$(awk -vS=$sample '($1==S){print $2}' SAMPLES.txt)
  ln -s /NextGenSeqData/project-data/daria/published_data/mES/Asifa_KANSL3_mES/bam/$sample.bam /NextGenSeqData/project-data/daria/published_data/mES/Asifa_KANSL3_mES/bam/links/$name.bam
done


mkdir -p ~/data/published_data/mES/Asifa_KANSL3_mES/macs2
screen -R macs
cd ~/data/published_data/mES/Asifa_KANSL3_mES/

for pair in KANSL3_mES_rep1-input_mES_rep1 KANSL3_mES_rep2-input_mES_rep2 KANSL3_NPC_rep1-input_NPC_rep1 KANSL3_NPC_rep2-input_NPC_rep2 ; do
   exp=$(echo $pair | awk '{split($1,a,"-");print a[1]}')
   input=$(echo $pair | awk '{split($1,a,"-");print a[2]}')
   echo -en "$exp\t$input"
   macs2 callpeak -t bam/links/$exp.bam -c bam/links/$input.bam -f BAM -n ${exp}_${input} --outdir macs2 --call-summits -g hs --keep-dup 1
done

# calling peaks on merged replicates
mkdir -p ~/data/published_data/mES/Asifa_KANSL3_mES/macs2
screen -R macs
cd ~/data/published_data/mES/Asifa_KANSL3_mES/

macs2 callpeak -t bam/links/KANSL3_mES_rep1.bam bam/links/KANSL3_mES_rep2.bam -c bam/links/input_mES_rep1.bam bam/links/input_mES_rep2.bam -f BAM -n KANSL3_mES_merged --outdir macs2 --call-summits -g hs --keep-dup 1 &

macs2 callpeak -t bam/links/KANSL3_NPC_rep1.bam bam/links/KANSL3_NPC_rep2.bam -c bam/links/input_NPC_rep1.bam bam/links/input_NPC_rep2.bam -f BAM -n KANSL3_NPC_merged --outdir macs2 --call-summits -g hs --keep-dup 1 &

######
# overlapping peaks from 2 reps and thresholding to get the final set

#cd ~/data/published_data/mES/Asifa_KANSL3_mES/macs2
#awk '($1!="#" && $1!="" && NR>1)' KANSL3_mES_rep1_input_mES_rep1_peaks.xls | wc -l # 14165
#wk '($1!="#" && $1!="" && NR>1)' KANSL3_mES_rep2_input_mES_rep2_peaks.xls | wc -l # 12704

# only of rep1 overlaps with rep2, no matter how large is the overlap
# 15585 loj
# threshold 6: 3332, still some peaks so-so
#bedtools intersect -loj \
#-a <(awk '($1!="#" && $1!="" && NR>1)' KANSL3_mES_rep1_input_mES_rep1_peaks.xls) \
#-b <(awk '($1!="#" && $1!="" && NR>1)' KANSL3_mES_rep2_input_mES_rep2_peaks.xls) | awk '($8>=6 && $18>=6)' 
# still gettting strange empty peaks: chr19 4078222 4078579 with high enrichment
# better to threshold on fold enrichments for both replicates


# to get a final peak based on both replicates (w/o merging them) decided to use only summits, extend to 200 bp and threshold (>=7) before overlaping 
# ideally: merge replicates and re-call peaks

#bedtools intersect -loj \
#-a <(awk -vOFS="\t" '($5>=7){print $1,$2-100,$3+100,$5}' KANSL3_mES_rep1_input_mES_rep1_summits.bed) \
#-b <(awk -vOFS="\t" '($5>=7){print $1,$2-100,$3+100,$5}' KANSL3_mES_rep2_input_mES_rep2_summits.bed) |  \
#awk '($5!="."){if($2<$6){summit=$2+100+int(($6-$2)/2)}else{summit=$6+100+int(($2-$6)/2)};print $1,summit,summit+1}' > KANSL3_mES_rep12.summits.txt

#bedtools intersect -loj \
#-a <(awk -vOFS="\t" '($5>=7){print $1,$2-100,$3+100,$5}' KANSL3_NPC_rep1_input_NPC_rep1_summits.bed) \
#-b <(awk -vOFS="\t" '($5>=7){print $1,$2-100,$3+100,$5}' KANSL3_NPC_rep2_input_NPC_rep2_summits.bed) |  \
#awk '($5!="."){if($2<$6){summit=$2+100+int(($6-$2)/2)}else{summit=$6+100+int(($2-$6)/2)};print $1,summit,summit+1}' > KANSL3_NPC_rep12.summits.txt

# overlap bw NPC and mES 

#bedtools intersect -u \
#-a <(awk -vOFS="\t" '{print $1,$2-250,$3+250}' KANSL3_mES_rep12.summits.txt) \
#-b <(awk -vOFS="\t" '{print $1,$2-250,$3+250}' KANSL3_NPC_rep12.summits.txt) | wc -l #2466 

####### numbers of overlaps for replicates with the same windows
cd ~/data/published_data/mES/Asifa_KANSL3_mES/macs2

bedtools intersect -u \
-a <(awk -vOFS="\t" '($5>=7){print $1,$2-250,$3+250}' KANSL3_mES_rep1_input_mES_rep1_summits.bed) \
-b <(awk -vOFS="\t" '($5>=7){print $1,$2-250,$3+250}' KANSL3_mES_rep2_input_mES_rep2_summits.bed) |  wc -l 

bedtools intersect -u \
-a <(awk -vOFS="\t" '($5>=7){print $1,$2-250,$3+250}' KANSL3_NPC_rep1_input_NPC_rep1_summits.bed) \
-b <(awk -vOFS="\t" '($5>=7){print $1,$2-250,$3+250}' KANSL3_NPC_rep2_input_NPC_rep2_summits.bed) |  wc -l 

####### overlaps mES and NPCs

# 6245; mES 9461; NPC 13156
# with $5>20: 2138 from mES 2592 and NPC 4841
bedtools intersect -u \
-a <(awk -vOFS="\t" '($5>=7){print $1,$2-250,$3+250}' KANSL3_mES_merged_summits.bed) \
-b <(awk -vOFS="\t" '($5>=7){print $1,$2-250,$3+250}' KANSL3_NPC_merged_summits.bed) |  wc -l 

# overlap between top peaks
bedtools intersect -u \
-a <(sort -k5,5gr KANSL3_mES_merged_summits.bed|  awk -vOFS="\t" '(NR<500){print $1,$2-250,$3+250}' ) \
-b <(sort -k5,5gr KANSL3_NPC_merged_summits.bed|  awk -vOFS="\t" '(NR<500){print $1,$2-250,$3+250}' ) |  wc -l 


######################## motif analysis with uniprobe hg19

# will have to overlap $peaks with all mapped motifs
#ls ~/data/data_genomes/motifs/uniprobe/mast_p03/hg19/NAME.txt.gz
# will count only once even if 2 overlaps
# no sorting since all files are sorted

screen -R motifs
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/
peaks=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
folderMotifs=~/data/data_genomes/motifs/uniprobe/mast_p03/hg19

touch motifs/counts.differential.txt
for motif in `ls -1 $folderMotifs `; do
  name=$(basename $motif .txt.gz)
  echo $name
  number=$(bedtools intersect -u -a $peaks -b <(zcat $folderMotifs/$motif) | wc -l)
  echo -en "$name\t$number\n" >> motifs/counts.differential.txt
done

# do the same for negatives used for dreme de novo


screen -R motifs
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/

# take the same number of negatives as positives (seems like A does it this way)
shuf -n 156 dreme/negative.txt > dreme/negative156.txt
negative=dreme/negative156.txt
folderMotifs=~/data/data_genomes/motifs/uniprobe/mast_p03/hg19

rm motifs/counts.negative.txt 
touch motifs/counts.negative.txt

# negatives are expanded to be 500 bp in length
# take the same number of negatives as positives (seems like A does it this way)
# good idea to do the same for positives - expand to a certain length [not done yet]
for motif in `ls -1 $folderMotifs `; do
  name=$(basename $motif .txt.gz)
  echo $name
  number=$(bedtools intersect -u -a <(awk -vOFS="\t" '($3-$2<=500){print $1,$2-int((500-($3-$2))/2),$3+int((500-($3-$2))/2),"NAME","1"}' $negative | sort -k1,1 -k2,2n) -b <(zcat $folderMotifs/$motif) | wc -l)
  echo -en "$name\t$number\n" >> motifs/counts.negative.txt
done

# calculating enrichments and p-values

totalDiff=$(awk 'END{print NR}' macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed)
totalNeg=$(awk 'END{print NR}' dreme/negative156.txt)

awk -vOFS="\t" -vFile="motifs/counts.differential.txt" -vD=$totalDiff -vN=$totalNeg 'BEGIN{while(getline<File){diffCounts[$1]=$2}}($1 in diffCounts){print $1,diffCounts[$1],D,$2,N}' motifs/counts.negative.txt | hyper.R -i - > motifs/counts.hyper.txt

# conclusion: 
# regadless the number of negatives, nothing is significant
# maybe too few motifs (180) and not representative 


#### HOCOMOCO v11 core set
 

screen -R motifs
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/
mkdir -p motifs/HOCOMOCOv11

peaks=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
negative=dreme/negative.txt
folderMotifs=~/data/data_genomes/motifs/HOCOMOCOv11/mast_p03/hg19

rm motifs/HOCOMOCOv11/counts.differential.txt 
touch motifs/HOCOMOCOv11/counts.differential.txt
for motif in `ls -1 $folderMotifs `; do
  name=$(basename $motif .txt.gz)
  #echo $name
  number=$(bedtools intersect -u -a $peaks -b <(zcat $folderMotifs/$motif) | wc -l)
  echo -en "$name\t$number\n" >> motifs/HOCOMOCOv11/counts.differential.txt
done &

rm motifs/HOCOMOCOv11/counts.negative.txt 
touch motifs/HOCOMOCOv11/counts.negative.txt

for motif in `ls -1 $folderMotifs `; do
  name=$(basename $motif .txt.gz)
  #echo $name
  number=$(bedtools intersect -u -a <(awk -vOFS="\t" '($3-$2<=500){print $1,$2-int((500-($3-$2))/2),$3+int((500-($3-$2))/2),"NAME","1"}' $negative | sort -k1,1 -k2,2n) -b <(zcat $folderMotifs/$motif) | wc -l)
  echo -en "$name\t$number\n" >> motifs/HOCOMOCOv11/counts.negative.txt
done &

totalDiff=$(awk 'END{print NR}' macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed)
totalNeg=$(awk 'END{print NR}' dreme/negative.txt)

awk -vOFS="\t" -vFile="motifs/HOCOMOCOv11/counts.differential.txt" -vD=$totalDiff -vN=$totalNeg 'BEGIN{while(getline<File){diffCounts[$1]=$2}}($1 in diffCounts){print $1,diffCounts[$1],D,$2,N,log((diffCounts[$1]/D)/($2/N))/log(2)}' motifs/HOCOMOCOv11/counts.negative.txt

awk -vOFS="\t" -vFile="motifs/HOCOMOCOv11/counts.differential.txt" -vD=$totalDiff -vN=$totalNeg 'BEGIN{while(getline<File){diffCounts[$1]=$2}}($1 in diffCounts){print $1,diffCounts[$1],D,$2,N}' motifs/HOCOMOCOv11/counts.negative.txt | hyper.R -i - > /tmp/counts.hyper.txt

# repeat motif search with random negatives (not CAGE)
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017

bedtools random -n 1000 -l 500 -g /home/daria/data/data_genomes/hg19/chrom.sizes.good > dreme/random.bed

screen -R motifs
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/
mkdir -p motifs/HOCOMOCOv11

random=dreme/random.bed
folderMotifs=~/data/data_genomes/motifs/HOCOMOCOv11/mast_p03/hg19

rm motifs/HOCOMOCOv11/counts.random.txt 
touch motifs/HOCOMOCOv11/counts.random.txt
for motif in `ls -1 $folderMotifs `; do
  name=$(basename $motif .txt.gz)
  #echo $name
  number=$(bedtools intersect -u -a $random -b <(zcat $folderMotifs/$motif) | wc -l)
  echo -en "$name\t$number\n" >> motifs/HOCOMOCOv11/counts.random.txt
done &

### counting enrichment relative to random and negative

totalDiff=$(awk 'END{print NR}' macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed)
totalNeg=$(awk 'END{print NR}' dreme/random.bed)

awk -vOFS="\t" -vFile="motifs/HOCOMOCOv11/counts.differential.txt" -vD=$totalDiff -vN=$totalNeg 'BEGIN{while(getline<File){diffCounts[$1]=$2}}($1 in diffCounts){print $1,diffCounts[$1],D,$2,N,log((diffCounts[$1]/D)/($2/N))/log(2)}' motifs/HOCOMOCOv11/counts.random.txt

awk -vOFS="\t" -vFile="motifs/HOCOMOCOv11/counts.differential.txt" -vD=$totalDiff -vN=$totalNeg 'BEGIN{while(getline<File){diffCounts[$1]=$2}}($1 in diffCounts){print $1,diffCounts[$1],D,$2,N}' motifs/HOCOMOCOv11/counts.random.txt | hyper.R -i -  | awk '{if($6<$7){p=$6}else{p=$7};print $1,$2,$3,$4,$5,log(($2/$3)/($4/$5))/log(2),p}' > /tmp/counts.hyper.random.txt

# for negative
totalDiff=$(awk 'END{print NR}' macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed)
totalNeg=$(awk 'END{print NR}' dreme/negative.txt)

awk -vOFS="\t" -vFile="motifs/HOCOMOCOv11/counts.differential.txt" -vD=$totalDiff -vN=$totalNeg 'BEGIN{while(getline<File){diffCounts[$1]=$2}}($1 in diffCounts){print $1,diffCounts[$1],D,$2,N,log((diffCounts[$1]/D)/($2/N))/log(2)}' motifs/HOCOMOCOv11/counts.negative.txt

awk -vOFS="\t" -vFile="motifs/HOCOMOCOv11/counts.differential.txt" -vD=$totalDiff -vN=$totalNeg 'BEGIN{while(getline<File){diffCounts[$1]=$2}}($1 in diffCounts){print $1,diffCounts[$1],D,$2,N}' motifs/HOCOMOCOv11/counts.negative.txt | hyper.R -i -  | awk '{if($6<$7){p=$6}else{p=$7};print $1,$2,$3,$4,$5,log(($2/$3)/($4/$5))/log(2),p}' > /tmp/counts.hyper.negative.txt

# negative vs random
totalDiff=$(awk 'END{print NR}' dreme/negative.txt)
totalNeg=$(awk 'END{print NR}' dreme/random.bed)

awk -vOFS="\t" -vFile="motifs/HOCOMOCOv11/counts.negative.txt" -vD=$totalDiff -vN=$totalNeg 'BEGIN{while(getline<File){diffCounts[$1]=$2}}($1 in diffCounts){print $1,diffCounts[$1],D,$2,N,log((diffCounts[$1]/D)/($2/N))/log(2)}' motifs/HOCOMOCOv11/counts.random.txt

awk -vOFS="\t" -vFile="motifs/HOCOMOCOv11/counts.negative.txt" -vD=$totalDiff -vN=$totalNeg 'BEGIN{while(getline<File){diffCounts[$1]=$2}}($1 in diffCounts){print $1,diffCounts[$1],D,$2,N}' motifs/HOCOMOCOv11/counts.random.txt | hyper.R -i -  | awk '{if($6<$7){p=$6}else{p=$7};print $1,$2,$3,$4,$5,log(($2/$3)/($4/$5))/log(2),p}' > /tmp/counts.hyper.negVSrandom.txt

# table for heatmaps: subselect KANSL vs rand motifs based on p-val and enrichment, and negative values for the same motifs --> KANSL vs negatives

cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/
awk '($7<=0.01 && ($6>=log(1.5)/log(2) || $6<=-log(1.5)/log(2))){print $1,$6}' /tmp/counts.hyper.random.txt | awk -vFile="/tmp/counts.hyper.negVSrandom.txt" 'BEGIN{while((getline<File)>0){motif[$1]=$6}}{if($1 in motif){print $1,$2,motif[$1]}else{print $1,$2,"NA"}}'  > /tmp/all.motifs.txt

awk '($7<=0.01){print $1,$6}' /tmp/counts.hyper.negative.txt > /tmp/peaksVSneg.txt

# tables for mouse data: 
cd ~/data/published_data/mES/Asifa_KANSL3_mES/motifs_counts

awk '($7<=0.01 && ($6>=log(1.5)/log(2) || $6<=-log(1.5)/log(2))){print $1,$6}' mES_vs_negative.txt | awk -vFile="NPC_vs_negative.txt" 'BEGIN{while((getline<File)>0){motif[$1]=$6}}{if($1 in motif){print $1,$2,motif[$1]}else{print $1,$2,"NA"}}'  > /tmp/all.motifs.mm9.txt

# will "center" on mES data
awk '($7<=0.01 && ($6>=log(1.5)/log(2) || $6<=-log(1.5)/log(2))){print $1,$6}' mES_vs_rand.txt | awk -vFile="NPC_vs_rand.txt" 'BEGIN{while((getline<File)>0){motif[$1]=$6}}{if($1 in motif){print $1,$2,motif[$1]}else{print $1,$2,"NA"}}' | awk -vFile="neg_vs_rand.txt" 'BEGIN{while((getline<File)>0){motif[$1]=$6}}{if($1 in motif){print $1,$2,$3,motif[$1]}else{print $1,$2,$3,"NA"}}'  > /tmp/mES_neg_NPC.motifs.mm9.txt


R # heatmaps
# from here: http://sebastianraschka.com/Articles/heatmaps_in_r.html
if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
   }

my_palette <- colorRampPalette(c("green4","green","white","violet","purple"))(100)

motifs.random<-read.table("/tmp/all.motifs.txt")
rnames <- motifs.random[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(motifs.random[,2:ncol(motifs.random)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names
colnames(mat_data) <- c("KANSL3 vs random","CAGE vs random")

pdf ("~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/plots/KANSL3_CAGE_random.pdf")
  heatmap.2(mat_data,
  #cellnote = mat_data,  # same data set for cell labels
  main = "KANSL3 and CAGE vs random", # heat map title
  #notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  #breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  cexCol=1,
  labRow = FALSE,
  #Rowv = FALSE, 
  #Colv = FALSE,
  ) 
dev.off()

# heatmap2 takes only 2 columns and 2 rows minimum. Need to trick it by dubplicating matrix

motifs.neg<-read.table("/tmp/peaksVSneg.txt")
rnames <- motifs.neg[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(motifs.neg[,2:ncol(motifs.neg)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names
colnames(mat_data) <- c("KANSL3 vs CAGE") 
mat_data<-cbind (mat_data, mat_data) 


pdf ("~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/plots/KANSL3vsCAGE.pdf")
  heatmap.2(mat_data,
  #cellnote = mat_data,  # same data set for cell labels
  #main = "Correlation", # heat map title
  #notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  #breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  Colv="NA",             # turn off column clustering
  cexCol=1,
  labRow = FALSE)    
dev.off()

# heatmaps for mouse data
motifs.random<-read.table("/tmp/mES_neg_NPC.motifs.mm9.txt")
rnames <- motifs.random[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(motifs.random[,2:ncol(motifs.random)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names
colnames(mat_data) <- c("mES vs random","NPC vs random", "CAGE vs random")

pdf ("~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/plots/mouse_mES_NPC_random.pdf")
  heatmap.2(mat_data,
  #cellnote = mat_data,  # same data set for cell labels
  main = "mES, NPC, negative vs random", # heat map title
  #notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  #breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  cexCol=1,
  labRow = FALSE,
  #Rowv = FALSE, 
  #Colv = FALSE,
  ) 
dev.off()

#motifs.random<-read.table("/tmp/all.motifs.mm9.txt")
#rnames <- motifs.random[,1]                            # assign labels in column 1 to "rnames"
#mat_data <- data.matrix(motifs.random[,2:ncol(motifs.random)])  # transform column 2-5 into a matrix
#rownames(mat_data) <- rnames                  # assign row names
#colnames(mat_data) <- c("mES vs CAGE","NPC vs CAGE")
#mat_data<-rbind (mat_data, mat_data) 

#df ("~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/plots/mES_NPC_negative.pdf")
 # heatmap.2(mat_data,
  #cellnote = mat_data,  # same data set for cell labels
  #main = "mES and NPC vs CAGE", # heat map title
  #notecol="black",      # change font color of cell labels to black
  #density.info="none",  # turns off density plot inside color legend
  #trace="none",         # turns off trace lines inside the heat map
  #margins =c(12,9),     # widens margins around plot
  #col=my_palette,       # use on color palette defined earlier
  #breaks=col_breaks,    # enable color transition at specified limits
  #dendrogram="row",     # only draw a row dendrogram
  #cexCol=1,
  #labRow = FALSE,
  #Rowv = FALSE, 
  #Colv = FALSE,
  ) 
#dev.off()


q()

##### checking if the motif analysis working properly with a published dataset

mkdir ~/data/published_data/ReMap2018/
cd ~/data/published_data/ReMap2018/

wget http://tagc.univ-mrs.fr/remap/download/MACS_lifted_hg19/all_tf/erg/ReMap2_erg_allPeaks_hg19.bed.tar.gz
gunzip ReMap2_erg_allPeaks_hg19.bed.tar.gz

## testing motif analysis with peaks for K562 cells 

screen -R motifs
cd ~/data/published_data/ReMap2018/

peaks=all_tf/erg/ReMap2_erg_allPeaks_hg19.bed
folderMotifs=~/data/data_genomes/motifs/HOCOMOCOv11/mast_p03/hg19
folder=~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017

rm $folder/test.txt
touch $folder/test.txt
for motif in `ls -1 $folderMotifs `; do
  name=$(basename $motif .txt.gz)
  #echo $name
  number=$(bedtools intersect -u -a <(awk '($4=="GSE23730.erg.k562")' $peaks) -b <(zcat $folderMotifs/$motif) | wc -l)
  echo -en "$name\t$number\n" >> $folder/test.txt
done &
# test looked very weired, no enrichment for Erg...... 

##### motifs in mouse data

screen -R motifs
cd ~/data/published_data/mES/Asifa_KANSL3_mES/

peaks=macs2/KANSL3_mES_merged_summits.bed
folderMotifs=~/data/data_genomes/motifs/HOCOMOCOv11/mast_p03/mm9
mkdir -p ~/data/published_data/mES/Asifa_KANSL3_mES/motifs_counts
folder=~/data/published_data/mES/Asifa_KANSL3_mES/motifs_counts

f=$folder/counts.mES.txt
rm $f
touch $f
for motif in `ls -1 $folderMotifs `; do
  name=$(basename $motif .txt.gz)
  #echo $name
  number=$(bedtools intersect -u -a <(awk -vOFS="\t" '($5>=20){print $1,$2-250,$3+250}' $peaks) -b <(zcat $folderMotifs/$motif) | wc -l)
  echo -en "$name\t$number\n" >> $f
done

screen -R random
cd ~/data/published_data/mES/Asifa_KANSL3_mES/

peaks=macs2/KANSL3_NPC_merged_summits.bed
folderMotifs=~/data/data_genomes/motifs/HOCOMOCOv11/mast_p03/mm9
mkdir -p ~/data/published_data/mES/Asifa_KANSL3_mES/motifs_counts
folder=~/data/published_data/mES/Asifa_KANSL3_mES/motifs_counts

f=$folder/counts.NPC.txt
rm $f
touch $f
for motif in `ls -1 $folderMotifs `; do
  name=$(basename $motif .txt.gz)
  #echo $name
  number=$(bedtools intersect -u -a <(awk -vOFS="\t" '($5>=20){print $1,$2-250,$3+250}' $peaks) -b <(zcat $folderMotifs/$motif) | wc -l)
  echo -en "$name\t$number\n" >> $f
done &

# making random (#2000) and negative controls 
cd ~/data/published_data/mES/Asifa_KANSL3_mES/
bedtools random -n 2000 -l 500 -g /home/daria/data/data_genomes/mm9/chrom.sizes.good > motifs_counts/random.bed

screen -R random
cd ~/data/published_data/mES/Asifa_KANSL3_mES/

peaks=motifs_counts/random.bed
folderMotifs=~/data/data_genomes/motifs/HOCOMOCOv11/mast_p03/mm9
mkdir -p ~/data/published_data/mES/Asifa_KANSL3_mES/motifs_counts
folder=~/data/published_data/mES/Asifa_KANSL3_mES/motifs_counts

f=$folder/counts.random.txt
rm $f
touch $f
for motif in `ls -1 $folderMotifs `; do
  name=$(basename $motif .txt.gz)
  #echo $name
  number=$(bedtools intersect -u -a <(sort -k1,1 -k2,2n $peaks) -b <(zcat $folderMotifs/$motif) | wc -l)
  echo -en "$name\t$number\n" >> $f
done &

# 2000 negative from TSSs: will take random p1 promoters, no data on expression levels in NPC and mES
# extend to 500 for overlapping

zcat ~/data/published_data/FANTOM5/TSS_mouse.bed.gz | awk -vOFS="\t" '{split($4,a,"@");win=250-int(($3-$2)/2);if(a[1]=="p1" && $1!="chrM"){print $1,$2-win,$3+win}}' | shuf -n 2000 | sort -k1,1 -k2,2n > motifs_counts/negative.txt

screen -R random
cd ~/data/published_data/mES/Asifa_KANSL3_mES/

peaks=motifs_counts/negative.txt
folderMotifs=~/data/data_genomes/motifs/HOCOMOCOv11/mast_p03/mm9
mkdir -p ~/data/published_data/mES/Asifa_KANSL3_mES/motifs_counts
folder=~/data/published_data/mES/Asifa_KANSL3_mES/motifs_counts

f=$folder/counts.negative.txt
rm $f
touch $f
for motif in `ls -1 $folderMotifs `; do
  name=$(basename $motif .txt.gz)
  #echo $name
  number=$(bedtools intersect -u -a <(sort -k1,1 -k2,2n $peaks) -b <(zcat $folderMotifs/$motif) | wc -l)
  echo -en "$name\t$number\n" >> $f
done &

#### comparisons
cd ~/data/published_data/mES/Asifa_KANSL3_mES
# negative vs random
totalDiff=$(awk '($5>=20)' macs2/KANSL3_mES_merged_summits.bed | wc -l)
totalDiffNPC=$(awk '($5>=20)' macs2/KANSL3_NPC_merged_summits.bed | wc -l)
totalNeg=$(awk 'END{print NR}' motifs_counts/negative.txt)
totalRand=$(awk 'END{print NR}' motifs_counts/random.bed)

file1=motifs_counts/counts.mES.txt
file2=motifs_counts/counts.NPC.txt
file3=motifs_counts/counts.negative.txt
file4=motifs_counts/counts.random.txt

#awk -vOFS="\t" -vFile=$file1 -vD=$totalDiff -vN=$totalNeg 'BEGIN{while(getline<File){diffCounts[$1]=$2}}(NR>2 && $1 in diffCounts){print $1,diffCounts[$1],D,$2,N,log((diffCounts[$1]/D)/($2/N))/log(2)}' $file3

awk -vOFS="\t" -vFile=$file1 -vD=$totalDiff -vN=$totalNeg 'BEGIN{while(getline<File){diffCounts[$1]=$2}}(NR>2 && $1 in diffCounts){print $1,diffCounts[$1],D,$2,N}' $file3 | hyper.R -i -  | awk '{if($6<$7){p=$6}else{p=$7};print $1,$2,$3,$4,$5,log(($2/$3)/($4/$5))/log(2),p}' > motifs_counts/mES_vs_negative.txt

awk -vOFS="\t" -vFile=$file1 -vD=$totalDiff -vN=$totalRand 'BEGIN{while(getline<File){diffCounts[$1]=$2}}(NR>2 && $1 in diffCounts){print $1,diffCounts[$1],D,$2,N}' $file4 | hyper.R -i -  | awk '{if($6<$7){p=$6}else{p=$7};print $1,$2,$3,$4,$5,log(($2/$3)/($4/$5))/log(2),p}' > motifs_counts/mES_vs_rand.txt

awk -vOFS="\t" -vFile=$file2 -vD=$totalDiffNPC -vN=$totalNeg 'BEGIN{while(getline<File){diffCounts[$1]=$2}}(NR>2 && $1 in diffCounts){print $1,diffCounts[$1],D,$2,N}' $file3 | hyper.R -i -  | awk '{if($6<$7){p=$6}else{p=$7};print $1,$2,$3,$4,$5,log(($2/$3)/($4/$5))/log(2),p}' > motifs_counts/NPC_vs_negative.txt

awk -vOFS="\t" -vFile=$file2 -vD=$totalDiffNPC -vN=$totalRand 'BEGIN{while(getline<File){diffCounts[$1]=$2}}(NR>2 && $1 in diffCounts){print $1,diffCounts[$1],D,$2,N}' $file4 | hyper.R -i -  | awk '{if($6<$7){p=$6}else{p=$7};print $1,$2,$3,$4,$5,log(($2/$3)/($4/$5))/log(2),p}' > motifs_counts/NPC_vs_rand.txt

awk -vOFS="\t" -vFile=$file3 -vD=$totalNeg -vN=$totalRand 'BEGIN{while(getline<File){diffCounts[$1]=$2}}(NR>2 && $1 in diffCounts){print $1,diffCounts[$1],D,$2,N}' $file4 | hyper.R -i -  | awk '{if($6<$7){p=$6}else{p=$7};print $1,$2,$3,$4,$5,log(($2/$3)/($4/$5))/log(2),p}' > motifs_counts/neg_vs_rand.txt

awk -vOFS="\t" -vFile=$file1 -vD=$totalDiff -vN=$totalDiffNPC 'BEGIN{while(getline<File){diffCounts[$1]=$2}}(NR>2 && $1 in diffCounts){print $1,diffCounts[$1],D,$2,N}' $file2 | hyper.R -i -  | awk '{if($6<$7){p=$6}else{p=$7};print $1,$2,$3,$4,$5,log(($2/$3)/($4/$5))/log(2),p}' > motifs_counts/mES_vs_NPC.txt

###################### H3K4me3 ChIP in WT and KD THP1

# copying all files from Desktop to the server

cd FASTQ_Generation*
for f in * ; do
echo -en "$f\n"
scp $f/ARD*.fastq.gz daria@tycho.sund.root.ku.dk:~/data/projects/BRIC_data/Alex/ChIP-seq/H3K4me3_THP1/raw/.
done

# on the server
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/H3K4me3_THP1
myTemp=~/data/temp_NOBACKUP
for f in `ls raw/*_R1_001.fastq.gz`; do   
   name=$(echo -en "${f%%.*}" | awk '{split($1,a,"_S");split(a[1],n,"ARD-"); print n[2]}')
   label=$(basename $name)
   echo -en "$label\n"
   zcat $f >> $myTemp/$label.merged.fastq
done

# number of reads: 
myTemp=~/data/temp_NOBACKUP
for name in K3SD2-H3K4me3-1 K3SD2-H3K4me3-2 THP1-H3K4me3-1 THP1-H3K4me3-2 THP1-input; do
   echo -en "$name\n"
   cat $myTemp/$name.merged.fastq | wc -l
done

#K3SD2-H3K4me3-1 
#99107360.00 24776840
#K3SD2-H3K4me3-2 
#94385728.00 23596432
#THP1-H3K4me3-1  
#103265540.00  25816385
#THP1-H3K4me3-2  
#118144520.00  29536130
#THP1-input  
#113926784.00  28481696

screen -R mapping
# do not have gz files - cat on files 
# will map to hg19 as all Fantom data is for hg19 + motifs
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/H3K4me3_THP1
myTemp=~/data/temp_NOBACKUP
hg19=/k/genomes/hg19/index/bowtie_canonical/hg19 # somehow index in my directory was not working
for name in K3SD2-H3K4me3-1 K3SD2-H3K4me3-2 THP1-H3K4me3-1 THP1-H3K4me3-2 THP1-input; do
  echo -en $name
  cat $myTemp/$name.merged.fastq | bowtie -p 5 -q -m 1 -v 3 --sam --best --strata --quiet $hg19 - > $myTemp/$name.sam
done

### making bam files
myTemp=~/data/temp_NOBACKUP
for sample in K3SD2-H3K4me3-1 K3SD2-H3K4me3-2 THP1-H3K4me3-1 THP1-H3K4me3-2 THP1-input; do
  echo $sample
  samtools view -Sb $myTemp/${sample}.sam > $myTemp/${sample}_nonSorted.bam
  samtools sort $myTemp/${sample}_nonSorted.bam > $myTemp/${sample}.bam
  samtools rmdup -s $myTemp/${sample}.bam $myTemp/${sample}.nodup.bam 
  samtools index $myTemp/${sample}.nodup.bam
  rm $myTemp/${sample}.sam $myTemp/${sample}_nonSorted.bam $myTemp/${sample}.bam
done

mkdir -p ~/data/projects/BRIC_data/Alex/ChIP-seq/H3K4me3_THP1/bam
newFolder=~/data/projects/BRIC_data/Alex/ChIP-seq/H3K4me3_THP1/bam
myTemp=~/data/temp_NOBACKUP
for sample in K3SD2-H3K4me3-1 K3SD2-H3K4me3-2 THP1-H3K4me3-1 THP1-H3K4me3-2 THP1-input; do
  mv $myTemp/${sample}.nodup.bam  $newFolder/.
  mv $myTemp/${sample}.nodup.bam.bai  $newFolder/.
done 

############# making bw files with macs2: OBS! macs2 does not normalize to the number of reads but rather to a control
# it would be needed when calling peaks, but not to represent data.... 
# better use bedtools genomecov and use a scale factor 1000000/mapped reads

#screen -R tracks
#cd ~/data/projects/BRIC_data/Alex/ChIP-seq/H3K4me3_THP1
#mkdir -p macs2; mkdir -p bw

# bam files without duplicates!
# something strange with chrM --> had to modify the entry for chrM to extend its length
#CHROMSIZE=chrom.sizes 
#for sample in K3SD2-H3K4me3-1 K3SD2-H3K4me3-2 THP1-H3K4me3-1 THP1-H3K4me3-2 THP1-input; do
 # echo -en "$sample\n"
#macs2 filterdup -i bam/${sample}.nodup.bam -o macs2/${sample}.bed
 # macs2 pileup -i macs2/${sample}.bed -o macs2/${sample}.pileup.bdg --extsize 140
  #bedGraphToBigWig macs2/${sample}.pileup.bdg $CHROMSIZE bw/${sample}.pileup.bw
#done

# making soft links and gettinf links for UCSC
#cd ~/web/projects/Alex
#folder=/NextGenSeqData/project-data/daria/projects/BRIC_data/Alex/ChIP-seq/H3K4me3_THP1/bw
#for file in $folder/*bw; do
 # ln -s $file .
#done

# making links for all files: 
#for file in ~/web/projects/Alex/*bw; do
 # name=$(basename $file)
 # echo -en "https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/$name\n"
#done

# with bedtools
# test
#file=K3SD2-H3K4me3-1.nodup
#hg19=~/data/data_genomes/hg19/chrom.sizes

#scaleFactor=$(samtools idxstats bam/$file.bam | awk '{total+=$3} END{print 1000000/total}')
#bedtools genomecov -ibam bam/$file.bam -bg -scale $scaleFactor -g $hg19 > bw/$file.bg
#-clip - If set just issue warning messages rather than dying if wig
 #                 file contains items off end of chromosome.
#wigToBigWig bw/$file.bg $hg19 bw/$file.normalised.bw

#### with bamCoverage tool --> will go with this one for now as only one line. Visually everything looked very simmilar
#### bw files are not so sharp, can play around with --extendReads and --smoothLength

cd ~/data/projects/BRIC_data/Alex/ChIP-seq/H3K4me3_THP1
file=K3SD2-H3K4me3-1.nodup
hg19=~/data/data_genomes/hg19/chrom.sizes # complained about chrM (is missing from good.sizes file)

for file in K3SD2-H3K4me3-1 K3SD2-H3K4me3-2 THP1-H3K4me3-1 THP1-H3K4me3-2 THP1-input; do
  scaleFactor=$(samtools idxstats bam/$file.nodup.bam | awk '{total+=$3} END{print 1000000/total}')
  #echo -en "$file\t$scaleFactor\n"
  # --centerReads shifts the center to make the enrichment site sharper 
  bamCoverage --bam bam/$file.nodup.bam -o bw/$file.SeqDepthNorm.bw --scaleFactor $scaleFactor --centerReads
done

#K3SD2-H3K4me3-1 0.0639786
#K3SD2-H3K4me3-2 0.0596021
#THP1-H3K4me3-1  0.0631394
#THP1-H3K4me3-2  0.0466284
#THP1-input  0.0439778

#K3SD2-H3K4me3-1 15630233
#K3SD2-H3K4me3-2 16777926
#THP1-H3K4me3-1  15837963
#THP1-H3K4me3-2  21446166
#THP1-input  22738739

cd ~/web/projects/Alex
folder=/NextGenSeqData/project-data/daria/projects/BRIC_data/Alex/ChIP-seq/H3K4me3_THP1/bw
# making soft links and gettinf links for UCSC
for file in $folder/*SeqDepthNorm.bw; do
  ln -s $file .
done

# making links for all files: 
for file in ~/web/projects/Alex/*SeqDepthNorm.bw; do
  name=$(basename $file)
  echo -en "https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/$name\n"
done

############################ new ChIPs 6 March 2018
# on the server 
mkdir ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018
mkdir ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018/raw

# copying files from Desktop to 
cd ~/Desktop/ARD-38966930/FASTQ*
for f in ARDseq* ; do
  echo -en "$f\n"
  scp $f/*.fastq.gz daria@tycho.sund.root.ku.dk:~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018/raw/.
done

cd ~/Desktop/Unind*/FASTQ*
for f in Undetermined* ; do
  echo -en "$f\n"
  scp $f/*.fastq.gz daria@tycho.sund.root.ku.dk:~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018/raw/.
done

# on the server
# merging lanes and making nice names for files
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018
mkdir ~/data/temp_NOBACKUP/newChIP # will remove afterwards
myTemp=~/data/temp_NOBACKUP/newChIP
for f in `ls raw/*_R1_001.fastq.gz`; do   
   name=$(echo -en "${f%%.*}" | awk '{split($1,a,"_S"); print a[1]}')
   label=$(basename $name)
   echo -en "$label\n"
   zcat $f >> $myTemp/$label.merged.fastq
done

# getting all the names 
for f in `ls raw/*_R1_001.fastq.gz`; do   
   name=$(echo -en "${f%%.*}" | awk '{split($1,a,"_S"); print a[1]}')
   label=$(basename $name)
   echo -en "$label\n" 
done > /tmp/names.txt
awk 

# number of reads: 
myTemp=~/data/temp_NOBACKUP/newChIP
for name in TK5a-1 K3K5a-1 K8K5a-1 TK5a-2 K3K5a-2 K8K5a-2 TK5u K3K5u K8K5u Ti K3i K8i TK16a-1 K3K16a-1 K8K16a-1 TK16a-2 K3K16a-2 K8K16a-2 K8K3me3; do
   echo -en "$name\n"
   cat $myTemp/$name.merged.fastq | wc -l
done

# mapping with bowtie
screen -R mapping
# do not have gz files - cat on files 
# will map to hg19 as all Fantom data is for hg19 + motifs
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018
myTemp=~/data/temp_NOBACKUP/newChIP
hg19=/k/genomes/hg19/index/bowtie_canonical/hg19 # somehow index in my directory was not working
for name in TK5a-1 K3K5a-1 K8K5a-1 TK5a-2 K3K5a-2 K8K5a-2 TK5u K3K5u K8K5u Ti; do
  echo -en $name
  cat $myTemp/$name.merged.fastq | bowtie -p 5 -q -m 1 -v 3 --sam --best --strata --quiet $hg19 - > $myTemp/$name.sam
done &

# the second batch
for name in K3i K8i TK16a-1 K3K16a-1 K8K16a-1 TK16a-2 K3K16a-2 K8K16a-2 K8K3me3; do
  echo -en $name
  cat $myTemp/$name.merged.fastq | bowtie -p 5 -q -m 1 -v 3 --sam --best --strata --quiet $hg19 - > $myTemp/$name.sam
done &

### making bam files
# in two batches
myTemp=~/data/temp_NOBACKUP/newChIP
for sample in TK5a-1 K3K5a-1 K8K5a-1 TK5a-2 K3K5a-2 K8K5a-2 TK5u K3K5u K8K5u Ti; do
  echo $sample
  samtools view -Sb $myTemp/${sample}.sam > $myTemp/${sample}_nonSorted.bam
  samtools sort $myTemp/${sample}_nonSorted.bam > $myTemp/${sample}.bam
  samtools rmdup -s $myTemp/${sample}.bam $myTemp/${sample}.nodup.bam 
  samtools index $myTemp/${sample}.nodup.bam
  rm $myTemp/${sample}.sam $myTemp/${sample}_nonSorted.bam $myTemp/${sample}.bam
done&

for sample in K3i K8i TK16a-1 K3K16a-1 K8K16a-1 TK16a-2 K3K16a-2 K8K16a-2 K8K3me3; do
  echo $sample
  samtools view -Sb $myTemp/${sample}.sam > $myTemp/${sample}_nonSorted.bam
  samtools sort $myTemp/${sample}_nonSorted.bam > $myTemp/${sample}.bam
  samtools rmdup -s $myTemp/${sample}.bam $myTemp/${sample}.nodup.bam 
  samtools index $myTemp/${sample}.nodup.bam
  rm $myTemp/${sample}.sam $myTemp/${sample}_nonSorted.bam $myTemp/${sample}.bam
done&


mkdir -p ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018/bam
newFolder=~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018/bam
myTemp=~/data/temp_NOBACKUP/newChIP
for sample in TK5a-1 K3K5a-1 K8K5a-1 TK5a-2 K3K5a-2 K8K5a-2 TK5u K3K5u K8K5u Ti K3i K8i TK16a-1 K3K16a-1 K8K16a-1 TK16a-2 K3K16a-2 K8K16a-2 K8K3me3; do
  mv $myTemp/${sample}.nodup.bam  $newFolder/.
  mv $myTemp/${sample}.nodup.bam.bai  $newFolder/.
done 

# remove fastq files and all content of the temp folder
rm $myTemp/*

# counting mapped reads
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018
for sample in TK5a-1 K3K5a-1 K8K5a-1 TK5a-2 K3K5a-2 K8K5a-2 TK5u K3K5u K8K5u Ti K3i K8i TK16a-1 K3K16a-1 K8K16a-1 TK16a-2 K3K16a-2 K8K16a-2 K8K3me3; do
  echo -en "$sample\t"
  samtools idxstats bam/$sample.nodup.bam | awk '{total+=$3}END{print total}'
done 

#### making bw files for visual inspection
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018
mkdir -p bw

hg19=~/data/data_genomes/hg19/chrom.sizes # complained about chrM (is missing from good.sizes file)
for file in TK5a-1 K3K5a-1 K8K5a-1 TK5a-2 K3K5a-2 K8K5a-2 TK5u K3K5u K8K5u Ti K3i K8i TK16a-1 K3K16a-1 K8K16a-1 TK16a-2 K3K16a-2 K8K16a-2 K8K3me3 ; do
  scaleFactor=$(samtools idxstats bam/$file.nodup.bam | awk '{total+=$3} END{print 1000000/total}')
  echo -en "$file\t$scaleFactor\n"
  # --centerReads shifts the center to make the enrichment site sharper 
  bamCoverage --bam bam/$file.nodup.bam -o bw/$file.SeqDepthNorm.bw --scaleFactor $scaleFactor --centerReads --extendReads 400 --binSize 10
done

# making links to tracks
cd ~/web/projects/Alex
folder=/NextGenSeqData/project-data/daria/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018/bw
# making soft links and gettinf links for UCSC
for file in $folder/*SeqDepthNorm.bw; do
  ln -s $file .
done

# making links for all files: 
for file in ~/web/projects/Alex/*SeqDepthNorm.bw; do
  name=$(basename $file)
  echo -en "https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/$name\n"
done

############################# average profiles over KANSL3 targets vs active CAGE vs inactive CAGE

# pull bam files for replicates and make bw files with similar settings: need for average profiles and for box plots
# only for H4K5 and H3K4me3
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018

for name in K3K5a TK5a K8K5a; do
  echo -en "$name\n"
  myTemp=~/data/temp_NOBACKUP/newChIP
  samtools merge $myTemp/${name}.bam bam/${name}-1.nodup.bam bam/${name}-2.nodup.bam
  samtools index $myTemp/${name}.bam
  scaleFactor=$(samtools idxstats $myTemp/${name}.bam | awk '{total+=$3} END{print 1000000/total}')
  bamCoverage --bam $myTemp/${name}.bam -o bw/${name}.SeqDepthNorm.merged.bw --scaleFactor $scaleFactor --centerReads --extendReads 400 --binSize 10
done

# the same for old data 
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/H3K4me3_THP1

for name in K3SD2-H3K4me3 THP1-H3K4me3; do
  echo -en "$name\n"
  myTemp=~/data/temp_NOBACKUP/newChIP
  samtools merge $myTemp/${name}.bam bam/${name}-1.nodup.bam bam/${name}-2.nodup.bam
  samtools index $myTemp/${name}.bam
  scaleFactor=$(samtools idxstats $myTemp/${name}.bam | awk '{total+=$3} END{print 1000000/total}')
  bamCoverage --bam $myTemp/${name}.bam -o bw/${name}.SeqDepthNorm.merged.bw --scaleFactor $scaleFactor --centerReads --extendReads 400 --binSize 10
done

# clean temdirectory 
rm ~/data/temp_NOBACKUP/newChIP/*bai ~/data/temp_NOBACKUP/newChIP/*bam 

# getting averages over bed files 
# still need to know how 0 are taken into consideration

cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018

# getting KANLS3 peaks, negative regions and random regions --> extend to 500 all of them (like peaks called by macs2)
myTemp=~/data/temp_NOBACKUP/newChIP/

awk -vOFS="\t" '(NR>1){summit=$2+int(($3-$2)/2);print $1,summit-250,summit+250,"KANSL3"NR}' ../KANSL3_Dec2017/macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed > $myTemp/KANSL3.bed
awk -vOFS="\t" '{summit=$2+int(($3-$2)/2);print $1,summit-250,summit+250,"neg"NR}' ../KANSL3_Dec2017/dreme/negative.txt> $myTemp/negative.bed
awk -vOFS="\t" '{summit=$2+int(($3-$2)/2);print $1,summit-250,summit+250,"rand"NR}' ../KANSL3_Dec2017/dreme/random.bed > $myTemp/rand.bed

for file in K3K5a-1.SeqDepthNorm.bw K3K5a-2.SeqDepthNorm.bw K3K5a.SeqDepthNorm.merged.bw K3K5u.SeqDepthNorm.bw TK5a-1.SeqDepthNorm.bw TK5a-2.SeqDepthNorm.bw TK5u.SeqDepthNorm.bw TK5a.SeqDepthNorm.merged.bw; do
  label=$(echo $file | awk '{split($1,a,"."); print a[1]}' )
  echo $label
  bigWigAverageOverBed bw/$file $myTemp/KANSL3.bed -bedOut=$myTemp/$label.KANSL3peaks.bed $myTemp/$label.KANSL3peaks.tab
  bigWigAverageOverBed bw/$file $myTemp/negative.bed -bedOut=$myTemp/$label.neg.bed $myTemp/$label.neg.tab
  bigWigAverageOverBed bw/$file $myTemp/rand.bed -bedOut=$myTemp/$label.rand.bed $myTemp/$label.rand.tab 
done

# for H3K4me3
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/H3K4me3_THP1
myTemp=~/data/temp_NOBACKUP/newChIP/
for file in THP1-H3K4me3-1.SeqDepthNorm.bw THP1-H3K4me3-2.SeqDepthNorm.bw THP1-H3K4me3.SeqDepthNorm.merged.bw K3SD2-H3K4me3-1.SeqDepthNorm.bw K3SD2-H3K4me3-2.SeqDepthNorm.bw K3SD2-H3K4me3.SeqDepthNorm.merged.bw; do
  label=$(echo $file | awk '{split($1,a,"."); print a[1]}' )
  echo $label
  bigWigAverageOverBed bw/$file $myTemp/KANSL3.bed -bedOut=$myTemp/$label.KANSL3peaks.bed $myTemp/$label.KANSL3peaks.tab
  bigWigAverageOverBed bw/$file $myTemp/negative.bed -bedOut=$myTemp/$label.neg.bed $myTemp/$label.neg.tab
  bigWigAverageOverBed bw/$file $myTemp/rand.bed -bedOut=$myTemp/$label.rand.bed $myTemp/$label.rand.tab 
done

# merging files for peaks, negative regions and random regions
# HAS TO RE-RUN before plotting
rm $myTemp/*forR.txt
for exp in KANSL3peaks neg rand; do
 touch $myTemp/$exp.forR.txt
 for file in TK5a-1.SeqDepthNorm.bw TK5a-2.SeqDepthNorm.bw TK5u.SeqDepthNorm.bw TK5a.SeqDepthNorm.merged.bw K3K5a-1.SeqDepthNorm.bw K3K5a-2.SeqDepthNorm.bw K3K5u.SeqDepthNorm.bw K3K5a.SeqDepthNorm.merged.bw; do
  label=$(echo $file | awk '{split($1,a,"."); print a[1]}' )
  echo $label
  paste $myTemp/$exp.forR.txt <(cut -f5 $myTemp/$label.$exp.bed) > $myTemp/temp.txt
  mv $myTemp/temp.txt $myTemp/$exp.forR.txt
done
done

# for H3K4me3 data --> re-wrote all K5ac :(
rm $myTemp/*forR.txt
for exp in KANSL3peaks neg rand; do
 touch $myTemp/$exp.forR.txt
 for file in THP1-H3K4me3-1.SeqDepthNorm.bw THP1-H3K4me3-2.SeqDepthNorm.bw THP1-H3K4me3.SeqDepthNorm.merged.bw K3SD2-H3K4me3-1.SeqDepthNorm.bw K3SD2-H3K4me3-2.SeqDepthNorm.bw K3SD2-H3K4me3.SeqDepthNorm.merged.bw ; do
  label=$(echo $file | awk '{split($1,a,"."); print a[1]}' )
  echo $label
  paste $myTemp/$exp.forR.txt <(cut -f5 $myTemp/$label.$exp.bed) > $myTemp/temp.txt
  mv $myTemp/temp.txt $myTemp/$exp.forR.txt
done
done

# final files
~/data/temp_NOBACKUP/newChIP/KANSL3peaks.forR.txt
~/data/temp_NOBACKUP/newChIP/rand.forR.txt 
~/data/temp_NOBACKUP/newChIP/neg.forR.txt

R 
source ("~/utils/Rutils/boxplot95.R")

file1<-read.table("~/data/temp_NOBACKUP/newChIP/KANSL3peaks.forR.txt")
file2<-read.table("~/data/temp_NOBACKUP/newChIP/rand.forR.txt")
file3<-read.table("~/data/temp_NOBACKUP/newChIP/neg.forR.txt")

#rownames(file1)<-c("TK5a-1","TK5a-2","TK5u","TK5a","K3K5a-1","K3K5a-2","K3K5u","K3K5a")

pdf("plots/boxplots_H4K5ac_WT_K3.pdf")
par(mfrow=c(4,1))

boxplot95(file1[,4],file2[,4],file3[,4],file1[,8],file2[,8],file3[,8],
  names=c("TK5a peaks","TK5a random","TK5a neg","K3K5a peaks","K3K5a random","K3K5a neg"),
  frame=F,main=c("merged replicates"),las=2)
boxplot95(file1[,3],file2[,3],file3[,3],file1[,7],file2[,7],file3[,7],
  names=c("TK5a peaks","TK5a random","TK5a neg","K3K5a peaks","K3K5a random","K3K5a neg"),
  frame=F,main=c("upstate ab"),las=2)
boxplot95(file1[,1],file2[,1],file3[,1],file1[,5],file2[,5],file3[,5],
  names=c("TK5a peaks","TK5a random","TK5a neg","K3K5a peaks","K3K5a random","K3K5a neg"),
  frame=F,main=c("rep 1"),las=2)
boxplot95(file1[,2],file2[,2],file3[,2],file1[,7],file2[,7],file3[,7],
  names=c("TK5a peaks","TK5a random","TK5a neg","K3K5a peaks","K3K5a random","K3K5a neg"),
  frame=F,main=c("rep 2"),las=2)
dev.off()

wilcox.test(file1[,8],file3[,8])
wilcox.test(file3[,4],file3[,8])

# for H3K4me3
source ("~/utils/Rutils/boxplot95.R")

file1<-read.table("~/data/temp_NOBACKUP/newChIP/KANSL3peaks.forR.txt")
file2<-read.table("~/data/temp_NOBACKUP/newChIP/rand.forR.txt")
file3<-read.table("~/data/temp_NOBACKUP/newChIP/neg.forR.txt")

# rows: WT-1 WT-2 WT K3-1 K3-2 K3

pdf("plots/boxplots_H3K4me3_WT_K3.pdf")
par(mfrow=c(3,1))

boxplot95(file1[,1],file2[,1],file3[,1],file1[,4],file2[,4],file3[,4],
  names=c("WT peaks","WT random","WT neg","K3 peaks","K3 random","K3 neg"),
  frame=F,main=c("rep 1, H3K4me3"),las=2)

boxplot95(file1[,2],file2[,2],file3[,2],file1[,5],file2[,5],file3[,5],
  names=c("WT peaks","WT random","WT neg","K3 peaks","K3 random","K3 neg"),
  frame=F,main=c("rep 2,H3K4me3"),las=2)

boxplot95(file1[,3],file2[,3],file3[,3],file1[,6],file2[,6],file3[,6],
  names=c("WT peaks","WT random","WT neg","K3 peaks","K3 random","K3 neg"),
  frame=F,main=c("merged, H3K4me3"),las=2)

dev.off()

q()

###### average profiles

# to be continued

################### OLD script from Alex's paper #########################
### possible script
#rm /tmp/out*
#names=~/data/projects/BRIC_data/Alex/CRISPR_CAGE_files/file_names.txt
#while read F  ; do
#  label=$(basename $F .bw)
#  echo $label
#  mark=~/data/published_data/ENCODE/human/K562/$F
#  bigWigAverageOverBed $mark /tmp/file.coord.txt -bedOut=/tmp/out.$label.bed /tmp/out.$label.tab
#done < $names

# average prfiles 

#bigWigToBedGraph ~/data/published_data/ENCODE/human/K562/wgEncodeSydhNsomeK562Sig.bw stdout > /tmp/mnase.bg &
#bigWigToBedGraph ~/data/published_data/ENCODE/human/K562/wgEncodeOpenChromDnaseK562Aln_2Reps.norm5.rawsignal.bw stdout > /tmp/dnase.bg &
#mnase=/tmp/mnase.bg
#dnase=/tmp/dnase.bg 
#win=1000
#awk -vWIN=$win '(NR>1 && $1=="chr1"){print $7,$8+int(($9-$8)/2)-WIN,$8+int(($9-$8)/2)+WIN,$6"_"$10}' $cage | sort -u | awk -vOFS="\t" '{chr[$4]=$1;start[$4]=$2;end[$4]=$3}END{for(k in chr){for(i=start[k];i<=end[k];i++)print chr[k],i,i+1,k}}' | sort -V -k1,1 -k2,2n > /tmp/cage.windows.bp.txt

#win=1000
#awk -vWIN=$win '(NR>1 && $6=="chr1"){print $6,$7+int(($8-$7)/2)-WIN,$7+int(($8-$7)/2)+WIN,$6"_"$7"_"$8}' $tss | sort -u | awk -vOFS="\t" '{chr[$4]=$1;start[$4]=$2;end[$4]=$3}END{for(k in chr){for(i=start[k];i<=end[k];i++)print chr[k],i,i+1,k}}' | sort -V -k1,1 -k2,2n > /tmp/tss.windows.bp.txt

#bedtools intersect -loj -a /tmp/cage.windows.bp.txt -b $mnase > /tmp/cage.windows.bp.mnase.txt &
#bedtools intersect -loj -a /tmp/cage.windows.bp.txt -b $dnase > /tmp/cage.windows.bp.dnase.txt &

#bedtools intersect -loj -a /tmp/tss.windows.bp.txt -b $mnase > /tmp/tss.windows.bp.mnase.txt &

#for data in mnase dnase; do
#  awk '{if($8=="."){print $1,$2,$3,$4,$5,$6,$7,"NA"}else{print $0}}' /tmp/cage.windows.bp.${data}.txt | awk -vOFS="\t" '{d=$4;density[d]=density[d]"\t"$8}END{for(summit in density){print summit,density[summit]}}' > /tmp/density.cage.${data}.txt
#done

#awk '{if($8=="."){print $1,$2,$3,$4,$5,$6,$7,"NA"}else{print $0}}' /tmp/tss.windows.bp.mnase.txt | awk -vOFS="\t" '{d=$4;density[d]=density[d]"\t"$8}END{for(summit in density){print summit,density[summit]}}' > /tmp/density.tss.mnase.txt
#R
#x<-read.table("/tmp/density.cage.mnase.txt",row.names=1)
#y<-read.table("/tmp/density.cage.dnase.txt",row.names=1)
#z<-read.table("/tmp/density.tss.mnase.txt",row.names=1)

#library("zoo")
#plot(rollapply(apply(x[,800:1500],2,mean,na.rm=TRUE),1,by=1,mean),type="l",col="red",frame.plot=F,ylab="Signal",xaxt="n",main="MNase",las=1,ylim=c(0,1))
#lines(rollapply(apply(z[,800:1500],2,mean,na.rm=TRUE),1,by=1,mean),type="l",col="blue",frame.plot=F,ylab="Signal",xaxt="n",main="MNase",las=1,ylim=c(0,1))
#plot(rollapply(apply(y,2,mean,na.rm=TRUE),1,by=1,mean),type="l",col="blue",frame.plot=F,las=1,ylim=c(0,40))


################ CAGE data combined with RNA-seq and KANSL3 peak status

# CAGE and KANSL3 distance distribution
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/

fileA=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
fileB=~/temp_NOBACKUP/FANTOM_table_THP1.rank.txt

# closest FANTOM5 peaks --> best solution, better than intersect
# temp file does not have a header! 
# will extend all peaks to 500 bp but CAGE peaks - will keep centers 

bedtools closest -a <(awk -vOFS="\t" '(NR>1){summit=$2+int(($3-$2)/2) ;print $1,summit-500,summit+500}' $fileA | sort -k1,1 -k2,2n) \
-b <(awk -vOFS="\t" '{summit=$2+int(($3-$2)/2); print $1,summit,summit+1}' $fileB | sort -k1,1 -k2,2n) -d | \
 awk '{print $7}' > /tmp/CAGE.distance.txt

# distance to 500 random peaks
bedtools random -n 500 -l 1 -g /home/daria/data/data_genomes/hg19/chrom.sizes.good | awk -vOFS="\t" '{print $1,$2,$3}' > /tmp/random.txt

bedtools closest -a <(awk -vOFS="\t" '{summit=$2+int(($3-$2)/2) ;print $1,summit-500,summit+500}' $fileA | sort -k1,1 -k2,2n /tmp/random.txt) \
-b <(awk -vOFS="\t" '{summit=$2+int(($3-$2)/2); print $1,summit,summit+1}' $fileB | sort -k1,1 -k2,2n) -d | \
 awk '{print $7}' > /tmp/rand.distance.txt

 # making a data frame for R
 # add pseudocaunt if 0
 cat <(awk '{if($1==0){count=1}else{count=$1};print log(count)/log(2),"KANSL3peaks"}' /tmp/CAGE.distance.txt) \
  <(awk '{if($1==0){count=1}else{count=$1};print log(count)/log(2),"random"}' /tmp/rand.distance.txt) > /tmp/data.txt

cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/
R
#library(ggplot2)
require(ggplot2)
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

cage<-read.table("/tmp/data.txt")


pdf("plots/distnance_500bp_centersCAGE.pdf")

ggplot(cage, aes(x = cage$V2, y = cage$V1, color=cage$V2)) +
  geom_boxplot() + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") +
  labs(title="Distance to CAGE-defined TSSs",x="", y = "log2 distance [bp]")

dev.off()

q()

########

cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/

fileA=macs2/differential/diff_c1_vs_c2_c3.0_cond2.bed
fileB=~/temp_NOBACKUP/FANTOM_table_THP1.rank.txt


# will overlap with all CAGE data and keep either the best ranked TSSs for each gene or p1=0
# 25699 of such CAGEs: p1=0 11310 --> 44% of all cage peaks are not expressed in THP1

# preselect CAGE peaks so they have either pRank1 or p1=0, but not both. Have to uniq genes
# give priority to pRank1 not to p1

awk '($10=="pRank1" || ($10==0 && $6=="p1")){print $0}' $fileB | \
awk '($6!="p1" && $10=="pRank1" || ($6=="p1" && $10==0)){print $5}' | uniq -c | awk '($1!=1){print $2}' > /tmp/exlude.txt
awk -vFile="/tmp/exlude.txt" 'BEGIN{while((getline<File)>0){name[$1]=$1}}{if($5 in name){if($10=="pRank1"){print $0}}else{if($10=="pRank1" || ($10==0 && $6=="p1")){print $0}}}' $fileB > ~/data/temp_NOBACKUP/cage_exluded.txt

# checking 
awk '($10=="pRank1" || ($10==0 && $6=="p1")){print $0}' ~/data/temp_NOBACKUP/cage_exluded.txt | \
awk '($6!="p1" && $10=="pRank1" || ($6=="p1" && $10==0)){print $5}' | uniq -c | sort -k1,1gr 

rm /tmp/exlude.txt

# redirect fileB
fileB=~/data/temp_NOBACKUP/cage_exluded.txt

# values for CAGE only from fresh sample, column 10 or #7 in the original file
bedtools closest -a <(awk -vOFS="\t" '(NR>1){summit=$2+int(($3-$2)/2) ;print $1,summit-500,summit+500}' $fileA | sort -k1,1 -k2,2n) \
-b <(awk -vOFS="\t" '{print $0}' $fileB | sort -k1,1 -k2,2n) -d | \
awk '($14==0){print $8,$10}' > ~/data/temp_NOBACKUP/kansl3_genes_cage.txt

# get 500 expressed CAGE peaks without KANSL3 peaks: negative control set
bedtools closest -a <(awk -vOFS="\t" '($10=="pRank1"){print $0}' $fileB | sort -k1,1 -k2,2n) \
-b <(awk -vOFS="\t" '(NR>1){summit=$2+int(($3-$2)/2) ;print $1,summit-500,summit+500}' $fileA | sort -k1,1 -k2,2n) -d | \
awk '($14>1000){print $5,$7}' | shuf -n 500  > ~/data/temp_NOBACKUP/genes_cage_noKansl.txt

# all CAGE data excluding 0
awk -vOFS="\t" '($10=="pRank1"){print $5,$7}' $fileB > ~/data/temp_NOBACKUP/ALL_cage_expressed.txt

# making a data frame for R

cat <(awk '{print $1,log($2)/log(2),"KANSL3target"}' ~/data/temp_NOBACKUP/kansl3_genes_cage.txt) \
<(awk '{print $1,log($2)/log(2),"notKANSL3target"}' ~/data/temp_NOBACKUP/genes_cage_noKansl.txt) \
<(awk '{print $1,log($2)/log(2),"expressedCAGE"}' ~/data/temp_NOBACKUP/ALL_cage_expressed.txt) > /tmp/dataframe.txt


R 
require(ggplot2)
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

data<-read.table("/tmp/dataframe.txt")


pdf("plots/CAGEvalues.pdf")

ggplot(data, aes(x = data$V3, y = data$V2, color=data$V3)) +
  geom_boxplot(outlier.shape = NA) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  labs(title="CAGE values",x="", y = "log2 CAGE normalized counts") +
  scale_color_manual(labels = c("All non-0 CAGE", "KANSL3 targets","random non-0 non-KANSL3 CAGE"), values = c("green", "red","blue"))

dev.off()

wilcox.test(subset(data[,2],data$V3 == "KANSL3target"),subset(data[,2],data$V3 == "expressedCAGE")) # very small
wilcox.test(subset(data[,2],data$V3 == "KANSL3target"),subset(data[,2],data$V3 == "notKANSL3target")) # 0.000164
wilcox.test(subset(data[,2],data$V3 == "expressedCAGE"),subset(data[,2],data$V3 == "notKANSL3target")) # 0.3039


q()

###### adding RNA-seq values to dataframe with KANSL3 targets and non-targets
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/KANSL3_Dec2017/

cat <(awk '{print $1,log($2)/log(2),"KANSL3target"}' ~/data/temp_NOBACKUP/kansl3_genes_cage.txt) \
<(awk '{print $1,log($2)/log(2),"notKANSL3target"}' ~/data/temp_NOBACKUP/genes_cage_noKansl.txt) > /tmp/dataframe.txt

rnaseq=../../RNA-seq/data_May2017/DESeq2/all_experiemnts.deseq2.merged.txt

awk -vFile=$rnaseq 'BEGIN{while((getline<File)>0){name[$2]=$7"\t"$8}}{if($1 in name){print $0,name[$1]}else{print $0,"NA","NA"}}' /tmp/dataframe.txt > /tmp/rnase.data.txt

# how many were recovered: 122 out of 160 CAGE peaks for KANSL3 targets; 393 out of 500 for notKANSL3targets
awk '($3=="KANSL3target" && $4!="NA")' /tmp/rnase.data.txt | wc -l 

# how many increasing by p-value? p=0.01 and p=0.05

for i in KANSL3target notKANSL3target; do
  for p in 0.01 0.05; do
    awk -vOFS="\t" -vP=$p -vI=$i '($4!="NA" && $5<=P && $3==I){if($4<0){down++}else{up++}}END{print I,"p="P,"down="down,"up="up}' /tmp/rnase.data.txt
  done
done

#KANSL3target  p=0.01  down=26 up=1 --> 21% down; 0.8% up
#KANSL3target  p=0.05  down=37 up=1 --> 30.3% down, 0.8% up
#notKANSL3target p=0.01  down=18 up=2 --> 4.5% down and 0.5% up
#notKANSL3target p=0.05  down=30 up=9 --> 7.6% down and 2.3% up

# rnaseq target that go down significantly: what are their CAGE values?
p=0.05
awk -vOFS="\t" -vP=$p '($4!="NA"){if($4<0 && $5<=P){print $0,"down"}else{print $0,"notUp"}}' /tmp/rnase.data.txt > /tmp/cage_down.txt

R 
require(ggplot2)
library(plyr)
# making data frame to make a cumulative barplots
df2 <- data.frame(target=rep(c("KANSL3 target genes", "non-target genes"), each=2),
                change=rep(c("up", "down"),2),
                percent=c(0.8, 30.3, 2.3, 7.6))
# sorting
df_sorted <- arrange(df2, target, rev(change)) 

# making cumulative
# + 0.2 to shift labels a bit upward for esthetics 
df_cumsum <- ddply(df_sorted, "target",
                   transform, label_ypos=cumsum(percent)+0.2)

pdf("plots/rnaseq_p0.05_percentage.pdf")
ggplot(data=df_cumsum, aes(x=target, y=percent, fill=change)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=percent), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()
dev.off()

# plotting boxplot and marking genes that go down
data<-read.table("/tmp/cage_down.txt")

# to make 5-95 boxplot
require(ggplot2)
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

pdf("plots/CAGE_values_down_targets_KANSL3.pdf")

sub<-subset(data,data$V6 == "down")
ggplot(data, aes(x = data$V3, y = data$V2, color=data$V3)) +
  geom_boxplot(outlier.shape = NA) + # plots sort of double plot with 'outliers'
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  geom_point(data=sub,aes(x = sub$V3, y = sub$V2, color=sub$V3),position = position_jitterdodge(dodge.width = 0.45)) +
  labs(title="CAGE values",x="", y = "log2 CAGE normalized counts") +
  scale_color_manual(labels = c("KANSL3 targets down", "non-KANSL3 targets down"), values = c("red", "blue"))

dev.off()

q()

####################################################
####################### use the scaling factor to nosmalise ChIP-seq tracks

#### making bw files for visual inspection
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018
mkdir -p bw

#> ScalingFactors(brd.tss) # Estimated using promoters
#     H4K16ac_WT_rep1 H4K16ac_KAT8_KD_rep1      H4K16ac_WT_rep2 
#           2.0155981            0.6233684            1.5584766 
#H4K16ac_KAT8_KD_rep2             Input_WT             Input_KD 
#           0.5160824            0.8836218            1.1198646 

hg19=~/data/data_genomes/hg19/chrom.sizes # complained about chrM (is missing from good.sizes file)
for name in Ti_0.8836218 K8i_1.1198646 TK16a-1_2.0155981 K8K16a-1_0.6233684 TK16a-2_1.5584766 K8K16a-2_0.5160824; do
  file=$(echo -en $name | awk '{split($1,a,"_"); print a[1]}')
  scaleFactor=$(echo -en $name | awk '{split($1,a,"_"); print a[2]}')
  echo -en "$file\t$scaleFactor\n"
  # --centerReads shifts the center to make the enrichment site sharper 
  bamCoverage --bam bam/$file.nodup.bam -o bw/$file.BRDnorm.bw --scaleFactor $scaleFactor --centerReads --extendReads 400 --binSize 10
done

# making links to tracks
cd ~/web/projects/Alex
folder=/NextGenSeqData/project-data/daria/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018/bw
# making soft links and gettinf links for UCSC
for file in $folder/*BRDnorm.bw; do
  ln -s $file .
done

# making links for all files: 
for file in ~/web/projects/Alex/*BRDnorm.bw; do
  name=$(basename $file)
  echo -en "https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/$name\n"
done

####################### making average profiles for K8 and K3 KD using deeptools

#> ScalingFactors(brd.tss) # Estimated using promoters --> USED FOR KAT8, WILL USE HERE AS WELL
#       H4K16ac_WT_rep1 H4K16ac_KANSL3_KD_rep1        H4K16ac_WT_rep2 
#             1.7743517              1.0239981              1.3677038 
#H4K16ac_KANSL3_KD_rep2               Input_WT               Input_KD 
#             0.7400619              0.6492857              0.8374623 

cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018
mkdir -p bw

hg19=~/data/data_genomes/hg19/chrom.sizes # complained about chrM (is missing from good.sizes file)
for name in Ti_0.6492857 K8i_0.8374623 TK16a-1_1.7743517 K3K16a-1_1.0239981 TK16a-2_1.3677038 K3K16a-2_0.7400619; do
  file=$(echo -en $name | awk '{split($1,a,"_"); print a[1]}')
  scaleFactor=$(echo -en $name | awk '{split($1,a,"_"); print a[2]}')
  echo -en "$file\t$scaleFactor\n"
  # --centerReads shifts the center to make the enrichment site sharper 
  bamCoverage --bam bam/$file.nodup.bam -o bw/$file.K3BRDnorm.bw --scaleFactor $scaleFactor --centerReads --extendReads 400 --binSize 10
done

# making links to tracks
cd ~/web/projects/Alex
folder=/NextGenSeqData/project-data/daria/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018/bw
# making soft links and gettinf links for UCSC
for file in $folder/*K3BRDnorm.bw; do
  ln -s $file .
done

# making links for all files: 
for file in ~/web/projects/Alex/*K3BRDnorm.bw; do
  name=$(basename $file)
  echo -en "https://bricweb.sund.ku.dk/bric-data/daria/projects/Alex/$name\n"
done

####### making average profiles with deeptools

cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018
mkdir -p average_profiles

# getting genes for hg19
cd ~/data/data_genomes/hg19
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

# the annotation file chromosome naming did not match genome fasta file -> added chr to the gtf file manually; remove non-canonical chromosomes

zcat Homo_sapiens.GRCh37.75.gtf.gz | awk '(NR<=5){print $0}' > /tmp/header.txt
zcat Homo_sapiens.GRCh37.75.gtf.gz | awk '(length($1)<=2){if($1!="MT"){print "chr"$0}else{start=substr($0,1,2);gsub("MT","chrM",start);end=substr($0,3); print start,end}}' > /tmp/file.txt
cat /tmp/header.txt /tmp/file.txt > ~/data/data_genomes/hg38/gtf/Homo_sapiens.GRCh37.75.canonical.gtf
rm /tmp/header.txt /tmp/file.txt 


# profiles over all genes in the genome
# bed file with all genes
GTF=~/data/data_genomes/hg38/gtf/Homo_sapiens.GRCh37.75.canonical.gtf
awk -vOFS="\t" '(NR>5 && $3=="gene"){gsub("\"","");gsub(";","");print $1,$4,$5,$10,"0",$7}' $GTF > ~/data/data_genomes/hg19/Homo_sapiens.GRCh37.75.canonical.bed

screen -R scale
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018
allGenes=~/data/data_genomes/hg19/Homo_sapiens.GRCh37.75.canonical.bed

computeMatrix scale-regions -R $allGenes -S bw/*.BRDnorm.bw -b 2000 -a 2000 \
  --regionBodyLength 2000 \
  --skipZeros -o average_profiles/matrix_K8_BRDnorm_gene1kb.gz \
  --outFileNameMatrix average_profiles/matrix_K8_BRDnorm_gene1kb.tab \
  --outFileSortedRegions average_profiles/regions_K8_BRDnorm_gene1kb.bed

plotHeatmap \
 -m average_profiles/matrix_K8_BRDnorm_gene1kb.gz \
 -out average_profiles/K8_BRDnorm_gene1kb.png \
 --colorMap YlGnBu \
 --regionsLabel 'all genes' \
 --heatmapHeight 15 \
 --plotTitle 'average profile' &

# for KANSL3

computeMatrix scale-regions -R $allGenes -S bw/*.K3BRDnorm.bw -b 2000 -a 2000 \
  --regionBodyLength 2000 \
  --skipZeros -o average_profiles/matrix_K3_BRDnorm_gene1kb.gz \
  --outFileNameMatrix average_profiles/matrix_K3_BRDnorm_gene1kb.tab \
  --outFileSortedRegions average_profiles/regions_K3_BRDnorm_gene1kb.bed

plotHeatmap \
 -m average_profiles/matrix_K3_BRDnorm_gene1kb.gz \
 -out average_profiles/K3_BRDnorm_gene1kb.png \
 --colorMap YlGnBu \
 --regionsLabel 'all genes' \
 --heatmapHeight 15 \
 --plotTitle 'average profile' &


# average profiles with pulished data in mouse
cd ~/data/data_genomes/mm9
wget ftp://ftp.ensembl.org/pub/release-65/gtf/mus_musculus/Mus_musculus.NCBIM37.65.gtf.gz

################################# BRD normalization for all samples at the same time

cd data

/usr/local/lib64/R-3.4.2/bin/R
library(Tightrope)
lst <- c(
"devtools", "stringr", "data.table", "triangle", "caTools", "ica", "FNN",
"igraph", "mixtools"
)
lapply(lst,function(x){library(x,character.only=TRUE)}) 
source("https://bioconductor.org/biocLite.R")
lst <- c(
"Rsamtools", "GenomicAlignments", "GenomicRanges", "GenomicFeatures",
"rtracklayer", "biomaRt"
)
lapply(lst,function(x){library(x,character.only=TRUE)})

DIR="~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018/bam/"

bam.files <- c(
paste(DIR,"TK16a-1.nodup.bam",sep=""),
paste(DIR,"K8K16a-1.nodup.bam",sep=""),
paste(DIR,"TK16a-2.nodup.bam",sep=""),
paste(DIR,"K8K16a-2.nodup.bam",sep=""),
paste(DIR,"K3K16a-1.nodup.bam",sep=""),
paste(DIR,"K3K16a-2.nodup.bam",sep=""),
paste(DIR,"Ti.nodup.bam",sep=""),
paste(DIR,"K8i.nodup.bam",sep="")
)

conditions <- c(
"H4K16ac_WT_rep1",
"H4K16ac_KAT8_KD_rep1",
"H4K16ac_WT_rep2",
"H4K16ac_KAT8_KD_rep2",
"H4K16ac_KANSL3_KD_rep1",
"H4K16ac_KANSL3_KD_rep2",
"Input_WT",
"Input_KD"
)

bam.flag <- scanBamFlag(isDuplicate = F, isUnmappedQuery = F)

chip <- c("H4K16ac_WT_rep1", "H4K16ac_KAT8_KD_rep1",
"H4K16ac_WT_rep2", "H4K16ac_KAT8_KD_rep2",
"H4K16ac_KANSL3_KD_rep1","H4K16ac_KANSL3_KD_rep2")
ctrl <- c("Input_WT", "Input_KD")

data("EGA75_human") # gene features
data("hg19_blacklist") 
data("CGI_hg19")

# making regions over which the counts will be made
TSS <- resize(EGA75_human$TSS, width = 4000, fix = "center", use.names = F)
# Define CGI regions as 4kb windows centered at CpG-islands annotated by UCSC
CGI <- resize(CGI_hg19, width = 4000, fix = "center", use.names = F)
CGI <- with(EGA75_human, CleanupGRanges(CGI, seqinfo, organism = organism$name))
# Make 4kb windows every 1kb over the whole genome
if(! LoadObj(WGW)) {
WGW <- GenomicTiling(EGA75_human$seqinfo, s = 1000, w = 4000)
SaveObj(WGW)
}

# Will re-run COUNTS for everything just in case

#if(! LoadObj(COUNTS)) { # Load previously saved COUNTS object if available
COUNTS <- new.env()
# Read counts over gene units
COUNTS$GNU <- ReadCountMatrix(
bam.files, EGA75_human$GNU, paired = F, names = conditions,
param = ScanBamParam(flag = bam.flag)
)
# Read counts over promoter regions
COUNTS$TSS <- ReadCountMatrix(
bam.files, TSS, paired = F, names = conditions,
param = ScanBamParam(flag = bam.flag)
)
# Read counts over CpG-islands
COUNTS$CGI <- ReadCountMatrix(
bam.files, CGI, paired = F, names = conditions,
param = ScanBamParam(flag = bam.flag)
)
# Read counts over the whole genome (4kb windows every 1kb)
COUNTS$WGW <- ReadCountMatrix(
bam.files, WGW, paired = F, names = conditions,
param = ScanBamParam(flag = bam.flag)
)
SaveObj(COUNTS) # Save COUNTS object as RData
#}

# will not build plots, because looked at the data before: only one cluster
# bdt parameter taken from the documentation by Benjamin

brd.gnu <- BRD(COUNTS$GNU, controls = ctrl, ncl=1, bdt = c(0.8, 0.05))
brd.tss <- BRD(COUNTS$TSS, controls = ctrl, ncl=1, bdt = c(0.8, 0.05))
brd.cgi <- BRD(COUNTS$CGI, controls = ctrl, ncl=1 ,bdt = c(0.8, 0.05))

ScalingFactors(brd.gnu) # Estimated using gene units
ScalingFactors(brd.tss) # Estimated using promoters
ScalingFactors(brd.cgi) # Estimated using CpG-islands

# conclusion: the scaling factors did not differ much: the variation is the same if I had run BRD several times
# will keep bw files like they were normalized before

q()

################################################################

cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_March_2018/average_profiles

# to see sample boundaries  look at matrix*.gz file, first line. Each row contains info for one gene and several groups
#"sample_boundaries":[0,600,1200,1800,2400,3000,3600]

# breaking matrix into 4 groups: "K8i.BRDnorm","K8K16a-1.BRDnorm","K8K16a-2.BRDnorm","Ti.BRDnorm","TK16a-1.BRDnorm","TK16a-2.BRDnorm"

R
# read.table can handle zipped files directly

x<-read.table("matrix_K8_BRDnorm_gene1kb.gz",skip=1)

# second way
library("zoo")
pdf("draft_K8_average_nolabels.pdf")
plot(rollapply(apply(x[,7:600],2,mean,na.rm=TRUE),15,by=4,median),type="l",col="red",frame.plot=F,ylab="Signal",xaxt="n",main="H4K16ac",las=1,ylim=c(0.3,1))
lines(rollapply(apply(x[,600:1200],2,mean,na.rm=TRUE),15,by=4,median),type="l",col="blue",frame.plot=F,ylab="Signal",xaxt="n",main="H4K16ac",las=1,ylim=c(0.3,1))
lines(rollapply(apply(x[,1200:1800],2,mean,na.rm=TRUE),15,by=4,median),type="l",col="deepskyblue4",frame.plot=F,ylab="Signal",xaxt="n",main="H4K16ac",las=1,ylim=c(0.3,1))
lines(rollapply(apply(x[,1800:2400],2,mean,na.rm=TRUE),15,by=4,median),type="l",col="firebrick4",frame.plot=F,ylab="Signal",xaxt="n",main="H4K16ac",las=1,ylim=c(0.3,1))
lines(rollapply(apply(x[,2400:3000],2,mean,na.rm=TRUE),15,by=4,median),type="l",col="chartreuse3",frame.plot=F,ylab="Signal",xaxt="n",main="H4K16ac",las=1,ylim=c(0.3,1))
lines(rollapply(apply(x[,3000:3600],2,mean,na.rm=TRUE),15,by=4,median),type="l",col="chartreuse4",frame.plot=F,ylab="Signal",xaxt="n",main="H4K16ac",las=1,ylim=c(0.3,1))
dev.off()


####################### data from August 2018: all ChIPs for chromatin marks re-done
# problems: some sequencing has very few reads

# manually transfer data to the created folder
mkdir -p ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_August_2018
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_August_2018
mkdir -p bam; mkdir -p bw; mkdir -p raw; mkdir -p raw2 # raw will contain data from 2 sequencing runs (10.08.18 and 09.09.18)

# manually transfering fastq files to the folder raw
# from the Desktop
# terminal locally
cd Desktop 
for f in `ls ARD090918_ChIP-seq-1/FASTQ_Generation_*/ARDseq*/*_R1_001.fastq.gz`; do
  scp -r $f daria@tycho.sund.root.ku.dk:~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_August_2018/raw/.
done

cd Desktop 
for f in `ls ARD100818_ChIP-seq-2/FASTQ_Generation_*/ARDseq*/*_R1_001.fastq.gz`; do
  scp -r $f daria@tycho.sund.root.ku.dk:~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_August_2018/raw2/.
done

# on the server
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_August_2018
mkdir -p ~/data/temp_NOBACKUP/August; myTemp=~/data/temp_NOBACKUP/August
for f in `ls raw/*_R1_001.fastq.gz`; do   
   name=$(echo -en "${f%%.*}" | awk '{split($1,a,"_S"); print a[1]}')
   label=$(basename $name)
   echo -en "$label\n"
   zcat $f >> $myTemp/$label.merged.fastq
done

# another sequencing run
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_August_2018
mkdir -p ~/data/temp_NOBACKUP/August2; myTemp=~/data/temp_NOBACKUP/August2
for f in `ls raw2/*_R1_001.fastq.gz`; do   
   name=$(echo -en "${f%%.*}" | awk '{split($1,a,"_S"); print a[1]}')
   label=$(basename $name)
   echo -en "$label\n"
   zcat $f >> $myTemp/$label.merged.fastq
done

# number of reads: 
#myTemp=~/data/temp_NOBACKUP/August
myTemp=~/data/temp_NOBACKUP/August2
#rm raw_reads_counts.txt
for f in `ls $myTemp/*merged.fastq`; do
  name=$(echo -en "${f%%.*}" | awk '{split($1,a,"_S"); print a[1]}')
  label=$(basename $name)
  count=$(cat $f | wc -l)
  echo -en "$label\t$count\n"
done >> raw_reads_counts.txt

screen -R mapping
# do not have gz files - cat on files 
# will map to hg19 as all Fantom data is for hg19 + motifs
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_August_2018/
myTemp=~/data/temp_NOBACKUP/August
hg19=/k/genomes/hg19/index/bowtie_canonical/hg19 # somehow index in my directory was not working
for f in `ls $myTemp/*merged.fastq`; do
  label=$(basename ${f%%.*})
  echo -en "$label\n"
  cat $f | bowtie -p 5 -q -m 1 -v 3 --sam --best --strata --quiet $hg19 - > $myTemp/$label.sam
done

myTemp=~/data/temp_NOBACKUP/August2
hg19=/k/genomes/hg19/index/bowtie_canonical/hg19 # somehow index in my directory was not working
for f in `ls $myTemp/*merged.fastq`; do
  label=$(basename ${f%%.*})
  echo -en "$label\n"
  cat $f | bowtie -p 5 -q -m 1 -v 3 --sam --best --strata --quiet $hg19 - > $myTemp/$label.sam
done

###### making bam files (duplicates removed)
### making bam files
screen -R bam
cd ~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_August_2018/
myTemp=~/data/temp_NOBACKUP/August
for f in `ls $myTemp/*sam`; do
  sample=$(basename ${f%%.*})
  echo -en "$sample\n"
  samtools view -Sb $myTemp/${sample}.sam > $myTemp/${sample}_nonSorted.bam
  samtools sort $myTemp/${sample}_nonSorted.bam > $myTemp/${sample}.bam
  samtools rmdup -s $myTemp/${sample}.bam $myTemp/${sample}.nodup.bam 
  samtools index $myTemp/${sample}.nodup.bam
  rm $myTemp/${sample}.sam $myTemp/${sample}_nonSorted.bam $myTemp/${sample}.bam
done
# moving everything to the bam folder
newFolder=~~/data/projects/BRIC_data/Alex/ChIP-seq/chromatin_marks_August_2018/bam
myTemp=~/data/temp_NOBACKUP/August
for f in `ls $myTemp/*sam`; do
  sample=$(basename ${f%%.*})
  echo -en "$sample\n"
  mv $myTemp/${sample}.nodup.bam  $newFolder/.
  mv $myTemp/${sample}.nodup.bam.bai  $newFolder/.
done 



























